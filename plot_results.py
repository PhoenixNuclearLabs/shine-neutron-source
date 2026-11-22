#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies. All rights reserved.
# ******************************************************************************

"""Plot the results of the Artemis models: the neutron and photon spectra in the cassette regions,
and the neutron and photon flux maps of the mesh tallies, from MCNP, OpenMC, or both.

    python plot_results.py            both
    python plot_results.py --mcnp     MCNP results (mcnp/artemis.im and mcnp/artemis.imsht)
    python plot_results.py --openmc   OpenMC results (the latest openmc/statepoint.*.h5)

The plots go in images/, with the geometry plots: artemis_spectrum_<particle>.png, with every
code given on the same axes, and artemis_<code>_<particle>_flux_<view>.png, for the pit and cavity
views of the plot command files mcnp/plot_*. Reading MCNP results needs mcnptools (which comes
with MCNP), and reading OpenMC results needs the openmc Python package.

Each flux figure has a panel for each plane of a view (xy and xz), all the same height, with a
shared color scale (FIGURES), the same for every code given: a vmax of None is the largest value
in the panels. Anything below vmin is shown in the vmin color instead of being left blank. There
is a contour line at each colorbar tick, labeled with its value and drawn on the colorbar too.
Material boundaries are drawn on top from the OpenMC geometry plot of the same window
(images/artemis_openmc_geom_<plane>_<view>.png, made by plot_geometry.py), which shows only
material boundaries (the MCNP plots also show the boundaries between cells of one material).
Each mesh bin is drawn as a solid block. The figures are about as wide as a GitHub page (at 150
dpi) and are saved as 256-color indexed PNGs, which look the same as full-color ones but are
much smaller.

All results are absolute: MCNP's are scaled by the total neutron yield with its fm cards, and
OpenMC's by the strength of the source. Fluxes are per cm^2 per second.
"""

import argparse
import glob
import io
import os
import sys
import time

import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.default'] = 'regular'

HERE = os.path.dirname(os.path.abspath(__file__))
MCNP_DIR = os.path.join(HERE, 'mcnp')
OPENMC_DIR = os.path.join(HERE, 'openmc')
IMAGES_DIR = os.path.join(HERE, 'images')

CASSETTE_CELLS = [106, 206]
SPECTRUM_TALLIES = {'neutron': 14, 'photon': 24}

# Mesh tallies: one for each plane and view, neutrons then photons
PLANES = ['xy', 'xz']
VIEWS = ['pit', 'cavity']
MESH_TALLIES = {(particle, plane, view): 34 + 10 * ((k * len(VIEWS) + m) * len(PLANES) + n)
                for k, particle in enumerate(['neutron', 'photon'])
                for m, view in enumerate(VIEWS) for n, plane in enumerate(PLANES)}

# Flux figures: (particle, view, title, units, color scale (type, vmin, vmax; None for the
# maximum), axis tick spacing [cm]). Log scales: the flux falls by orders of magnitude from the
# source in the chamber to the cassettes (cavity views) and through the water (pit views)
FIGURES = [
    ('neutron', 'pit'   , 'Neutron flux'           , 'n/cm$^2$-s'        , ('log', 1e6, None), 100),
    ('neutron', 'cavity', 'Neutron flux (cavities)', 'n/cm$^2$-s'        , ('log', 1e8, None), 20),
    ('photon' , 'pit'   , 'Photon flux'            , r'$\gamma$/cm$^2$-s', ('log', 1e6, None), 100),
    ('photon' , 'cavity', 'Photon flux (cavities)' , r'$\gamma$/cm$^2$-s', ('log', 1e8, None), 20),
]

# Contour line colors: the first N colors of CONTOUR_COLORS_2 (N = number of contour levels), in the
# order they appear in CONTOUR_COLORS_1
CONTOUR_COLORS_1 = ['w', 'gold', 'r', 'b', 'g', 'm', 'c', 'darkred', 'darkblue', 'darkgreen',
                    'purple', 'darkcyan', 'dodgerblue', 'yellow', 'lime', 'aqua', 'ivory']
CONTOUR_COLORS_2 = ['r', 'b', 'g', 'yellow', 'w', 'm', 'darkblue', 'aqua', 'gold', 'darkred', 'c',
                    'darkgreen', 'purple', 'lime', 'darkcyan', 'ivory', 'dodgerblue']


def contour_colors(n):
    """Colors for n contour levels."""
    return [x for x in CONTOUR_COLORS_1 if x in CONTOUR_COLORS_2[:n]]


# ------------------------------------------------------------------------------
# Reading results
# ------------------------------------------------------------------------------
class MCNPResults:
    """MCNP results: the mctal (spectra) and meshtal (mesh tallies) files of mcnp6 n=artemis.i."""

    name = 'MCNP'

    def __init__(self):
        import mcnptools
        self.mctal = mcnptools.Mctal(os.path.join(MCNP_DIR, 'artemis.im'))
        self.meshtal = mcnptools.Meshtal(os.path.join(MCNP_DIR, 'artemis.imsht'))

    def spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the cassettes."""
        t = self.mctal.GetTally(SPECTRUM_TALLIES[particle])
        upper = np.array(t.GetEBins())  # upper edges; the first is the top of a bin from 0
        flux = np.array([[t.GetValue(c, 0, 0, 0, 0, 0, e, 0) for e in range(len(upper))]
                         for c in range(len(CASSETTE_CELLS))]).mean(axis=0)
        return upper, flux[1:]

    def mesh(self, particle, plane, view):
        """Return (extent, flux [/cm^2-s], rows along the vertical axis) of a mesh tally."""
        t = self.meshtal.GetTally(MESH_TALLIES[particle, plane, view])
        edges = [np.array(t.GetXRBounds()), np.array(t.GetYZBounds()), np.array(t.GetZTBounds())]
        shape = [len(e) - 1 for e in edges]
        flux = np.array([[[t.GetValue(i, j, k, 0) for k in range(shape[2])]
                          for j in range(shape[1])] for i in range(shape[0])])
        return in_plane(edges, flux, plane)


class OpenMCResults:
    """OpenMC results: the latest statepoint in openmc/, with the volumes of the model's cells."""

    name = 'OpenMC'

    def __init__(self):
        import openmc
        files = sorted(glob.glob(os.path.join(OPENMC_DIR, 'statepoint.*.h5')),
                       key=lambda f: int(f.split('.')[-2]))
        if not files:
            raise SystemExit(f'No statepoint in {OPENMC_DIR}: run openmc/artemis.py first')
        self.statepoint = openmc.StatePoint(files[-1])
        # Cell volumes, from the model (the tallies hold track lengths, not fluxes)
        sys.path.insert(0, OPENMC_DIR)
        import artemis
        cells = artemis.get_geometry(artemis.get_materials()).get_all_cells()
        self.volumes = [cells[c].volume for c in CASSETTE_CELLS]

    def spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the cassettes."""
        import openmc
        t = self.statepoint.get_tally(id=SPECTRUM_TALLIES[particle])
        edges = t.find_filter(openmc.EnergyFilter).values * 1e-6
        track = t.mean.reshape(len(CASSETTE_CELLS), -1)
        return edges, (track / np.array(self.volumes)[:, None]).mean(axis=0)

    def mesh(self, particle, plane, view):
        """Return (extent, flux [/cm^2-s], rows along the vertical axis) of a mesh tally."""
        import openmc
        t = self.statepoint.get_tally(id=MESH_TALLIES[particle, plane, view])
        mesh = t.find_filter(openmc.MeshFilter).mesh
        nx, ny, nz = mesh.dimension
        edges = [np.linspace(lo, hi, n + 1)
                 for lo, hi, n in zip(mesh.lower_left, mesh.upper_right, mesh.dimension)]
        volume = np.prod((np.array(mesh.upper_right) - np.array(mesh.lower_left))
                         / np.array(mesh.dimension))
        # Mesh bins run with x fastest, then y, then z
        flux = t.mean.ravel().reshape(nz, ny, nx).transpose(2, 1, 0) / volume
        return in_plane(edges, flux, plane)


def in_plane(edges, flux, plane):
    """Reduce a mesh tally one bin thick across the plane (flux[x, y, z]) to the panel's extent
    and its flux, with rows along the plane's vertical axis."""
    h, v = 'xyz'.index(plane[0]), 'xyz'.index(plane[1])
    flux = np.take(flux, 0, axis=3 - h - v)
    return [edges[h][0], edges[h][-1], edges[v][0], edges[v][-1]], flux.T


# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------
def geometry_edges(plane, view, width=1):
    """Material boundaries from the OpenMC geometry plot of a window, as an RGBA image (black where
    neighboring pixels differ in color), or None if there is no plot. The edges are widened by width
    pixels on each side, so they stay dark when the image is shrunk to the panel size."""
    filename = os.path.join(IMAGES_DIR, f'artemis_openmc_geom_{plane}_{view}.png')
    if not os.path.exists(filename):
        return None
    img = plt.imread(filename)[:, :, :3]
    edges = np.zeros(img.shape[:2], dtype=bool)
    edges[:, 1:] |= np.any(img[:, 1:] != img[:, :-1], axis=2)
    edges[1:, :] |= np.any(img[1:, :] != img[:-1, :], axis=2)
    wide = edges.copy()
    for shift in range(1, width + 1):
        wide[:, shift:] |= edges[:, :-shift]
        wide[:, :-shift] |= edges[:, shift:]
        wide[shift:, :] |= edges[:-shift, :]
        wide[:-shift, :] |= edges[shift:, :]
    rgba = np.zeros(img.shape[:2] + (4,))
    rgba[wide, 3] = 1.0
    return rgba


def plot_flux(results):
    """One figure per code and entry of FIGURES, with a panel for each plane."""
    for particle, view, title, units, (scale, vmin, vmax), tick in FIGURES:
        maps = {r.name: [r.mesh(particle, plane, view) for plane in PLANES] for r in results}
        top = max(np.max(flux) for panels in maps.values() for _, flux in panels)
        norm = (LogNorm if scale == 'log' else Normalize)(vmin, top if vmax is None else vmax)
        for name, panels in maps.items():
            # The panels all have the same height, each with its own aspect ratio
            height = 3.0  # [in]
            widths = [height * (e[1] - e[0]) / (e[3] - e[2]) for e, _ in panels]
            figsize = (sum(widths) + 0.6 * len(panels) + 0.9, height + 1.3)
            fig, axs = plt.subplots(1, len(panels), figsize=figsize, layout='constrained',
                                    width_ratios=widths)
            for ax, plane, (extent, flux) in zip(axs, PLANES, panels):
                im = ax.imshow(np.maximum(flux, norm.vmin), origin='lower', extent=extent,
                               norm=norm, cmap='Spectral_r', interpolation='nearest')
                edges = geometry_edges(plane, view)
                if edges is not None:
                    ax.imshow(edges, origin='upper', extent=extent)
                ax.set_title(f'{plane} plane')
                ax.set_xlabel(f'{plane[0]} [cm]')
                ax.set_ylabel(f'{plane[1]} [cm]')
                ax.xaxis.set_major_locator(MultipleLocator(tick))
                ax.yaxis.set_major_locator(MultipleLocator(tick))
            cbar = fig.colorbar(im, ax=axs, label=units, shrink=0.8, pad=0.01)

            # Contour lines at the colorbar ticks (not including the ends of the scale), one
            # color per tick, the same in every panel, drawn on the colorbar too, and labeled
            # with their values
            levels = [x for x in cbar.get_ticks() if norm.vmin < x < norm.vmax]
            levels = levels[:len(CONTOUR_COLORS_2)]
            colors = contour_colors(len(levels))
            for ax, (extent, flux) in zip(axs, panels):
                cs = ax.contour(flux, levels=levels, colors=colors, linewidths=0.8, origin='lower',
                                extent=extent)
                ax.clabel(cs, fontsize=6, fmt=lambda x: f'{x:.3g}'.replace('e+', 'e'))
            cbar.add_lines(cs)
            fig.suptitle(f'Artemis: {title}, {name}')
            fig.get_layout_engine().set(h_pad=0.2, wspace=0)
            save(fig, f'artemis_{name.lower()}_{particle}_flux_{view}.png', indexed=True)


def plot_spectra(results):
    """The cassette spectra, per unit lethargy, with every code on the same axes."""
    for particle in ['neutron', 'photon']:
        fig, ax = plt.subplots(figsize=(8, 6), layout='constrained')
        for r in results:
            edges, flux = r.spectrum(particle)
            ax.stairs(flux / np.log(edges[1:] / edges[:-1]), edges, label=r.name)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(f'{particle.capitalize()} energy [MeV]')
        ax.set_ylabel(f'{particle.capitalize()} flux per unit lethargy '
                      f'[{"n" if particle == "neutron" else "$\\gamma$"}/cm$^2$-s]')
        ax.set_title(f'Artemis: {particle} spectrum in the cassette regions '
                     f'(mean of cells {CASSETTE_CELLS[0]} and {CASSETTE_CELLS[1]})')
        ax.grid(which='major', alpha=0.4)
        ax.legend()
        save(fig, f'artemis_spectrum_{particle}.png')


def save(fig, name, indexed=False):
    """Save a figure in images/, at 150 dpi; indexed: as a 256-color indexed PNG (about a quarter
    the size of a full-color one, and looks the same: the colormap has 256 colors)."""
    from PIL import Image
    filename = os.path.join(IMAGES_DIR, name)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150)
    plt.close(fig)
    image = Image.open(buf).convert('RGB')
    if indexed:
        image = image.quantize(256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE)
    # On Windows, another program (an editor's preview, antivirus) can hold the file for a moment
    for attempt in range(5):
        try:
            image.save(filename, optimize=True)
            break
        except OSError:
            if attempt == 4:
                raise
            time.sleep(1)
    print('Wrote', os.path.relpath(filename, HERE))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--mcnp', action='store_true', help='MCNP results only')
    parser.add_argument('--openmc', action='store_true', help='OpenMC results only')
    args = parser.parse_args()
    both = not (args.mcnp or args.openmc)

    results = []
    if args.mcnp or both:
        results.append(MCNPResults())
    if args.openmc or both:
        results.append(OpenMCResults())
    plot_spectra(results)
    plot_flux(results)


if __name__ == '__main__':
    main()
