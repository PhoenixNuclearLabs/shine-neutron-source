#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies. All rights reserved.
# ******************************************************************************

"""Plot the results of the Artemis models: the neutron and photon spectra in the cassette regions,
and the neutron and photon flux maps of the mesh tallies, from MCNP, OpenMC, or both.

    python plot_results.py            both
    python plot_results.py --mcnp     MCNP results (mcnp/artemis.im and mcnp/artemis.imsht)
    python plot_results.py --openmc   OpenMC results (the latest openmc/statepoint.*.h5)

The plots go in images/, with the geometry plots: artemis_spectrum_<particle>.png, the spectra in
the cassette regions and at their peak flux, with every code given on the same axes, and flux maps
for the pit and cavity views of the plot command files mcnp/plot_*: the neutron flux in total and
in four energy groups (NEUTRON_GROUPS), artemis_<code>_neutron_flux_total_<view>.png and
artemis_<code>_neutron_flux_<group>_<view>.png, and the photon flux,
artemis_<code>_photon_flux_<view>.png. The peak neutron flux in each
cassette region (PEAK_TALLIES), in each energy group and in total, is printed as a table. Reading
MCNP results needs mcnptools (which comes with MCNP), and reading OpenMC results needs the openmc
Python package.

Each flux figure has a panel for each plane of a view (xy and xz), all the same height, with a
shared color scale (FIGURES), the same for every code given: a vmax of None is the largest value
in the panels. Anything below vmin is shown in the vmin color instead of being left blank. There
is a contour line at each colorbar tick, labeled with its value and drawn on the colorbar too.
Material boundaries are drawn on top from the OpenMC geometry plot of the same window
(images/artemis_openmc_geom_<plane>_<view>.png, made by plot_geometry.py), which shows only
material boundaries (the MCNP plots also show the boundaries between cells of one material). The
cassette regions, which are air like the cavities around them, have their open sides dashed.
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
import warnings

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"

HERE = os.path.dirname(os.path.abspath(__file__))
MCNP_DIR = os.path.join(HERE, "mcnp")
OPENMC_DIR = os.path.join(HERE, "openmc")
IMAGES_DIR = os.path.join(HERE, "images")

CASSETTE_CELLS = [106, 206]
SPECTRUM_TALLIES = {"neutron": 14, "photon": 24}

# Mesh tallies: one for each plane and view, neutrons then photons
PLANES = ["xy", "xz"]
VIEWS = ["pit", "cavity"]
MESH_TALLIES = {
    (particle, plane, view): 34 + 10 * ((k * len(VIEWS) + m) * len(PLANES) + n)
    for k, particle in enumerate(["neutron", "photon"])
    for m, view in enumerate(VIEWS)
    for n, plane in enumerate(PLANES)
}
# Peak flux in each cassette region: one-bin mesh tallies over a plate against the cavity wall, at
# the height of the peak flux: neutrons in the four groups, for each side, and the neutron and
# photon spectra (same bins as the cassette spectra), -x then +x
PEAK_TALLIES = {"-x": 114, "+x": 124}
PEAK_SPECTRUM_TALLIES = {"neutron": [134, 144], "photon": [154, 164]}

# Energy groups of the neutron mesh tallies (emesh in the MCNP deck, an energy filter in the
# OpenMC model): name (in the file names) and title of each, in order of energy
NEUTRON_GROUPS = [
    ("thermal", "below 0.5 eV"),
    ("epithermal", "from 0.5 eV to 0.1 MeV"),
    ("fast", "from 0.1 to 10 MeV"),
    ("high", "above 10 MeV"),
]
GROUP = {name: g for g, (name, _) in enumerate(NEUTRON_GROUPS)}

# Flux figures: (particle, energy group (None for all energies), view, color scale (type, vmin,
# vmax; None for the maximum), axis tick spacing [cm]). Log scales: the flux falls by orders of
# magnitude from the source in the chamber to the cassettes (cavity views) and through the water
# (pit views)
# fmt: off
FIGURES = [
    ('neutron', None        , 'pit'   , ('log', 1e6, None), 100),
    ('neutron', None        , 'cavity', ('log', 1e8, None), 20),
    ('neutron', 'thermal'   , 'pit'   , ('log', 1e6, None), 100),
    ('neutron', 'thermal'   , 'cavity', ('log', 1e7, None), 20),
    ('neutron', 'epithermal', 'pit'   , ('log', 1e6, None), 100),
    ('neutron', 'epithermal', 'cavity', ('log', 1e7, None), 20),
    ('neutron', 'fast'      , 'pit'   , ('log', 1e6, None), 100),
    ('neutron', 'fast'      , 'cavity', ('log', 1e7, None), 20),
    ('neutron', 'high'      , 'pit'   , ('log', 1e6, None), 100),
    ('neutron', 'high'      , 'cavity', ('log', 1e7, None), 20),
    ('photon' , None        , 'pit'   , ('log', 1e6, None), 100),
    ('photon' , None        , 'cavity', ('log', 1e8, None), 20),
]
# fmt: on
UNITS = {"neutron": "n/cm$^2$-s", "photon": r"$\gamma$/cm$^2$-s"}

# Cassette region of the +x cavity (the -x one is its mirror image), as in the models [cm]: a box,
# outside the cavity wall around the target chamber. It is air, like the cavity around it, so the
# flux maps draw its open sides as dashed lines (its other sides are walls)
CASSETTE_X = (1.27, 51.181)
CASSETTE_Y = 20.32  # half-width
CASSETTE_Z = (16.51, 94.615)
CASSETTE_R = 7.06501  # cavity wall around the target chamber
CASSETTE_STYLE = {"color": "dimgray", "linestyle": (0, (4, 2)), "linewidth": 1.2}

# Contour levels on log scales, in each decade: 1 and 3 for the whole-geometry (pit) views, and 1, 2,
# and 5 for the closer cavity views (linear scales have them at the colorbar ticks)
CONTOUR_STEPS = {"pit": [1, 3], "cavity": [1, 2, 5]}

# SHINE brand colors for the spectrum plots: the blue family and the grays (the product-line
# colors are only for series that are the product lines)
SHINE = {
    "dark blue": "#0093FF",
    "shine blue": "#00CCFF",
    "dark gray": "#53565A",
    "mid gray": "#B1B3B3",
    "light gray": "#CCCCCC",
    "ink": "#001027",
}
SPECTRUM_COLORS = {
    "peak": {"MCNP": SHINE["dark blue"], "OpenMC": SHINE["shine blue"]},
    "cassette": {"MCNP": SHINE["dark gray"], "OpenMC": SHINE["mid gray"]},
}

# Contour line colors: the first N colors of CONTOUR_COLORS_2 (N = number of contour levels), in the
# order they appear in CONTOUR_COLORS_1
CONTOUR_COLORS_1 = [
    "w",
    "gold",
    "r",
    "b",
    "g",
    "m",
    "c",
    "darkred",
    "darkblue",
    "darkgreen",
    "purple",
    "darkcyan",
    "dodgerblue",
    "yellow",
    "lime",
    "aqua",
    "ivory",
]
CONTOUR_COLORS_2 = [
    "r",
    "b",
    "g",
    "yellow",
    "w",
    "m",
    "darkblue",
    "aqua",
    "gold",
    "darkred",
    "c",
    "darkgreen",
    "purple",
    "lime",
    "darkcyan",
    "ivory",
    "dodgerblue",
]


def power_of_ten(x):
    """A contour label in the colorbar's style: 2×10⁸ (10⁹ for 1×10⁹), in mathtext."""
    mantissa, exponent = f"{x:.2e}".split("e")
    mantissa, exponent = f"{float(mantissa):g}", int(exponent)
    return f"$10^{{{exponent}}}$" if mantissa == "1" else f"${mantissa}\\times10^{{{exponent}}}$"


def contour_colors(n):
    """Colors for n contour levels."""
    return [x for x in CONTOUR_COLORS_1 if x in CONTOUR_COLORS_2[:n]]


# ------------------------------------------------------------------------------
# Reading results
# ------------------------------------------------------------------------------
class MCNPResults:
    """MCNP results: the mctal (spectra) and meshtal (mesh tallies) files of mcnp6 n=artemis.i."""

    name = "MCNP"

    def __init__(self):
        import mcnptools

        self.mctal = mcnptools.Mctal(os.path.join(MCNP_DIR, "artemis.im"))
        self.meshtal = mcnptools.Meshtal(os.path.join(MCNP_DIR, "artemis.imsht"))

    def spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the cassettes."""
        t = self.mctal.GetTally(SPECTRUM_TALLIES[particle])
        upper = np.array(t.GetEBins())  # upper edges; the first is the top of a bin from 0
        flux = np.array(
            [
                [t.GetValue(c, 0, 0, 0, 0, 0, e, 0) for e in range(len(upper))]
                for c in range(len(CASSETTE_CELLS))
            ]
        ).mean(axis=0)
        return upper, flux[1:]

    def peak_spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the plates."""
        # The edges from the cassette spectrum, which has the same bins (the meshtal file rounds
        # them to 3 digits)
        upper = np.array(self.mctal.GetTally(SPECTRUM_TALLIES[particle]).GetEBins())
        flux = []
        for n in PEAK_SPECTRUM_TALLIES[particle]:
            t = self.meshtal.GetTally(n)
            assert len(t.GetEBounds()) == len(upper) + 1
            flux.append([t.GetValue(0, 0, 0, e) for e in range(len(upper))])
        return upper, np.mean(flux, axis=0)[1:]

    def peak(self, side):
        """Return (flux [/cm^2-s], relative error) in each energy group and in total (last) of a
        peak flux tally."""
        t = self.meshtal.GetTally(PEAK_TALLIES[side])
        bins = range(len(t.GetEBounds()))  # the groups, then the total
        return (np.array([t.GetValue(0, 0, 0, e) for e in bins]),
                np.array([t.GetError(0, 0, 0, e) for e in bins]))  # fmt: skip

    def mesh(self, particle, plane, view):
        """Return (edges along x, y, and z [cm], flux[x, y, z, energy group] [/cm^2-s]) of a mesh
        tally."""
        t = self.meshtal.GetTally(MESH_TALLIES[particle, plane, view])
        edges = [np.array(t.GetXRBounds()), np.array(t.GetYZBounds()), np.array(t.GetZTBounds())]
        shape = [len(e) - 1 for e in edges] + [len(t.GetEBounds()) - 1]
        flux = np.array(
            [
                [
                    [[t.GetValue(i, j, k, e) for e in range(shape[3])] for k in range(shape[2])]
                    for j in range(shape[1])
                ]
                for i in range(shape[0])
            ]
        )
        return edges, flux


class OpenMCResults:
    """OpenMC results: the latest statepoint in openmc/, with the volumes of the model's cells."""

    name = "OpenMC"

    def __init__(self):
        import openmc

        files = sorted(
            glob.glob(os.path.join(OPENMC_DIR, "statepoint.*.h5")),
            key=lambda f: int(f.split(".")[-2]),
        )
        if not files:
            raise SystemExit(f"No statepoint in {OPENMC_DIR}: run openmc/artemis.py first")
        self.statepoint = openmc.StatePoint(files[-1])
        # Cell volumes, from the model (the tallies hold track lengths, not fluxes)
        sys.path.insert(0, OPENMC_DIR)
        import artemis

        # The statepoint's summary has already made the model's materials, surfaces, and cells:
        # making them again warns that their IDs are in use
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", openmc.IDWarning)
            cells = artemis.get_geometry(artemis.get_materials()).get_all_cells()
        self.volumes = [cells[c].volume for c in CASSETTE_CELLS]

    def spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the cassettes."""
        import openmc

        t = self.statepoint.get_tally(id=SPECTRUM_TALLIES[particle])
        edges = t.find_filter(openmc.EnergyFilter).values * 1e-6
        track = t.mean.reshape(len(CASSETTE_CELLS), -1)
        return edges, (track / np.array(self.volumes)[:, None]).mean(axis=0)

    def peak_spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the plates."""
        import openmc

        flux = []
        for n in PEAK_SPECTRUM_TALLIES[particle]:
            t = self.statepoint.get_tally(id=n)
            mesh = t.find_filter(openmc.MeshFilter).mesh
            volume = np.prod(np.array(mesh.upper_right) - np.array(mesh.lower_left))
            flux.append(t.mean.ravel() / volume)
        edges = t.find_filter(openmc.EnergyFilter).values * 1e-6
        return edges, np.mean(flux, axis=0)

    def peak(self, side):
        """Return (flux [/cm^2-s], relative error) in each energy group and in total (last) of a
        peak flux tally. The total's error treats the groups as independent."""
        import openmc

        t = self.statepoint.get_tally(id=PEAK_TALLIES[side])
        mesh = t.find_filter(openmc.MeshFilter).mesh
        volume = np.prod(np.array(mesh.upper_right) - np.array(mesh.lower_left))
        mean, sd = t.mean.ravel() / volume, t.std_dev.ravel() / volume
        mean = np.append(mean, mean.sum())
        sd = np.append(sd, np.sqrt(np.sum(sd**2)))
        return mean, np.divide(sd, mean, out=np.zeros_like(sd), where=mean > 0)

    def mesh(self, particle, plane, view):
        """Return (edges along x, y, and z [cm], flux[x, y, z, energy group] [/cm^2-s]) of a mesh
        tally."""
        import openmc

        t = self.statepoint.get_tally(id=MESH_TALLIES[particle, plane, view])
        mesh = t.find_filter(openmc.MeshFilter).mesh
        nx, ny, nz = mesh.dimension
        edges = [
            np.linspace(lo, hi, n + 1)
            for lo, hi, n in zip(mesh.lower_left, mesh.upper_right, mesh.dimension)
        ]
        volume = np.prod(
            (np.array(mesh.upper_right) - np.array(mesh.lower_left)) / np.array(mesh.dimension)
        )
        # The filters are mesh, particle, and energy (if any), the last varying fastest; mesh bins
        # run with x fastest, then y, then z
        flux = t.mean.ravel().reshape(nz, ny, nx, -1).transpose(2, 1, 0, 3) / volume
        return edges, flux


def in_plane(edges, flux, plane, group):
    """Reduce a mesh tally one bin thick across the plane (flux[x, y, z, energy group]) to the
    panel's extent and its flux in one energy group (the sum of them all for None), with rows along
    the plane's vertical axis."""
    flux = flux.sum(axis=3) if group is None else flux[:, :, :, GROUP[group]]
    h, v = "xyz".index(plane[0]), "xyz".index(plane[1])
    flux = np.take(flux, 0, axis=3 - h - v)
    return [edges[h][0], edges[h][-1], edges[v][0], edges[v][-1]], flux.T


# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------
def geometry_edges(plane, view, width=1):
    """Material boundaries from the OpenMC geometry plot of a window, as an RGBA image (black where
    neighboring pixels differ in color), or None if there is no plot. The edges are widened by width
    pixels on each side, so they stay dark when the image is shrunk to the panel size."""
    filename = os.path.join(IMAGES_DIR, f"artemis_openmc_geom_{plane}_{view}.png")
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


def cassette_sides(plane, at):
    """The open sides of both cassette regions in a plane at at [cm] across it, as (horizontal,
    vertical) coordinate pairs of line segments."""
    (x0, x1), y, (z0, z1) = CASSETTE_X, CASSETTE_Y, CASSETTE_Z
    segments = []
    if plane == "xy" and z0 < at < z1:
        # The sides at +-y and the outer end; the inner end is the cavity walls
        segments = [([x0, x1], [y, y]), ([x0, x1], [-y, -y]), ([x1, x1], [-y, y])]
    elif plane == "xz" and abs(at) < y:
        # The bottom, top, and outer end; the inner end is the cavity wall around the chamber
        inner = np.sqrt(CASSETTE_R**2 - at**2) if abs(at) < CASSETTE_R else x0
        segments = [([inner, x1], [z0, z0]), ([inner, x1], [z1, z1]), ([x1, x1], [z0, z1])]
    return [([sign * x for x in xs], ys) for xs, ys in segments for sign in [1, -1]]


def plot_flux(results):
    """One figure per code and entry of FIGURES, with a panel for each plane."""
    tallies = {}  # each mesh tally, read once
    for particle, group, view, (scale, vmin, vmax), tick in FIGURES:
        maps = {}
        for r in results:
            maps[r.name] = []
            for plane in PLANES:
                key = r.name, particle, plane, view
                if key not in tallies:
                    tallies[key] = r.mesh(particle, plane, view)
                maps[r.name].append(in_plane(*tallies[key], plane, group))
        # Where each plane is, across it: the middle of the mesh's single bin across the plane
        across = {}
        for plane in PLANES:
            edges = tallies[results[0].name, particle, plane, view][0]
            across[plane] = np.mean(edges[3 - "xyz".index(plane[0]) - "xyz".index(plane[1])][:2])
        title = f"{particle.capitalize()} flux"
        if group is not None:
            title += " " + dict(NEUTRON_GROUPS)[group]
        if view == "cavity":
            title += " (cavities)"
        units = UNITS[particle]
        # The neutron flux in total is named like its energy groups; the photon flux is all total
        if group is not None:
            suffix = f"{group}_{view}"
        else:
            suffix = f"total_{view}" if particle == "neutron" else view
        top = max(np.max(flux) for panels in maps.values() for _, flux in panels)
        norm = (LogNorm if scale == "log" else Normalize)(vmin, top if vmax is None else vmax)
        for name, panels in maps.items():
            # The panels all have the same height, each with its own aspect ratio
            height = 3.0  # [in]
            widths = [height * (e[1] - e[0]) / (e[3] - e[2]) for e, _ in panels]
            figsize = (sum(widths) + 0.6 * len(panels) + 0.9, height + 1.3)
            fig, axs = plt.subplots(
                1, len(panels), figsize=figsize, layout="constrained", width_ratios=widths
            )
            for ax, plane, (extent, flux) in zip(axs, PLANES, panels):
                im = ax.imshow(
                    np.maximum(flux, norm.vmin),
                    origin="lower",
                    extent=extent,
                    norm=norm,
                    cmap="Spectral_r",
                    interpolation="nearest",
                )
                edges = geometry_edges(plane, view)
                if edges is not None:
                    ax.imshow(edges, origin="upper", extent=extent)
                for xs, ys in cassette_sides(plane, across[plane]):
                    ax.plot(xs, ys, **CASSETTE_STYLE)
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
                ax.set_title(f"{plane} plane")
                ax.set_xlabel(f"{plane[0]} [cm]")
                ax.set_ylabel(f"{plane[1]} [cm]")
                ax.xaxis.set_major_locator(MultipleLocator(tick))
                ax.yaxis.set_major_locator(MultipleLocator(tick))
            cbar = fig.colorbar(im, ax=axs, label=units, shrink=0.8, pad=0.01)

            # Contour lines at CONTOUR_STEPS (not including the ends of the scale), one color per
            # level, the same in every panel, drawn on the colorbar too, and labeled with their
            # values
            if scale == "log":
                decades = range(
                    int(np.floor(np.log10(norm.vmin))), int(np.ceil(np.log10(norm.vmax))) + 1
                )
                levels = [step * 10.0**d for d in decades for step in CONTOUR_STEPS[view]]
            else:
                levels = cbar.get_ticks()
            levels = [x for x in levels if norm.vmin < x < norm.vmax]
            levels = levels[: len(CONTOUR_COLORS_2)]
            colors = contour_colors(len(levels))
            for ax, (extent, flux) in zip(axs, panels if levels else []):
                cs = ax.contour(
                    flux,
                    levels=levels,
                    colors=colors,
                    linewidths=0.8,
                    origin="lower",
                    extent=extent,
                )
                ax.clabel(cs, fontsize=6, fmt=power_of_ten)
            if levels:
                cbar.add_lines(cs)
            fig.suptitle(f"Artemis: {title}, {name}")
            fig.get_layout_engine().set(h_pad=0.2, wspace=0)
            save(fig, f"artemis_{name.lower()}_{particle}_flux_{suffix}.png", indexed=True)


def plot_spectra(results):
    """The spectra in the cassette regions and at their peak flux (on the plates), per unit
    lethargy, with every code on the same axes, in SHINE colors: blues for the plates, grays for the
    cassette regions, the darker shade for MCNP."""
    for particle in ["neutron", "photon"]:
        fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
        for peak in [True, False]:  # the cassette means (smoother) on top
            for r in results:
                edges, flux = r.peak_spectrum(particle) if peak else r.spectrum(particle)
                label = (
                    f"{r.name}, peak (5 x 5 x 0.1 cm plates)"
                    if peak
                    else f"{r.name}, cassette regions (cells {CASSETTE_CELLS[0]} and "
                    f"{CASSETTE_CELLS[1]})"
                )
                ax.stairs(
                    # Bins with no score are left out (gaps), instead of drops to the axis
                    np.where(flux > 0, flux, np.nan) / np.log(edges[1:] / edges[:-1]),
                    edges,
                    label=label,
                    color=SPECTRUM_COLORS["peak" if peak else "cassette"][r.name],
                    linewidth=1.0,
                    baseline=None,  # no closing lines down to zero at the ends of each run
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"{particle.capitalize()} energy [MeV]", color=SHINE["dark gray"])
        ax.set_ylabel(
            f"{particle.capitalize()} flux per unit lethargy "
            f"[{'n' if particle == 'neutron' else '$\\gamma$'}/cm$^2$-s]",
            color=SHINE["dark gray"],
        )
        ax.set_title(
            f"Artemis: {particle} spectrum in the cassette regions (mean of -x and +x)",
            color=SHINE["ink"],
        )
        ax.tick_params(which="both", colors=SHINE["dark gray"])
        for spine in ax.spines.values():
            spine.set_color(SHINE["dark gray"])
        ax.grid(which="major", color=SHINE["light gray"], linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(loc="upper left")
        save(fig, f"artemis_spectrum_{particle}.png")


def print_peak(results):
    """A table of the peak neutron flux in each cassette region, by energy group and in total."""
    columns = [title for _, title in NEUTRON_GROUPS] + ["total"]
    print("\nPeak neutron flux in the cassette regions [n/cm^2-s] (relative error)")
    print(f"{'':14}" + "".join(f"{c:>24}" for c in columns))
    for r in results:
        for side in PEAK_TALLIES:
            flux, error = r.peak(side)
            cells = "".join(f"{f'{v:.3e} ({e:.3f})':>24}" for v, e in zip(flux, error))
            print(f"{f'{r.name} {side}':14}" + cells)


def save(fig, name, indexed=False):
    """Save a figure in images/, at 150 dpi; indexed: as a 256-color indexed PNG (about a quarter
    the size of a full-color one, and looks the same: the colormap has 256 colors)."""
    from PIL import Image

    filename = os.path.join(IMAGES_DIR, name)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    image = Image.open(buf).convert("RGB")
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
    print("Wrote", os.path.relpath(filename, HERE))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mcnp", action="store_true", help="MCNP results only")
    parser.add_argument("--openmc", action="store_true", help="OpenMC results only")
    args = parser.parse_args()
    both = not (args.mcnp or args.openmc)

    results = []
    if args.mcnp or both:
        results.append(MCNPResults())
    if args.openmc or both:
        results.append(OpenMCResults())
    plot_spectra(results)
    plot_flux(results)
    print_peak(results)


if __name__ == "__main__":
    main()
