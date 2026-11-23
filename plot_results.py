#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies, LLC. MIT License: see LICENSE.
# ******************************************************************************

"""Plot the results of the Artemis models: the neutron and photon spectra in the cassette regions,
and the neutron and photon flux maps of the mesh tallies, from MCNP, OpenMC, or both.

    python plot_results.py            both
    python plot_results.py --mcnp     MCNP results (mcnp/artemis.im and mcnp/artemis.imsht)
    python plot_results.py --openmc   OpenMC results (the latest openmc/statepoint.*.h5)

The plots go in images/, with the geometry plots: artemis_spectrum_<particle>.png, the spectra in
the cassette regions and at their peak flux, with every code given on the same axes, and maps for
the pit and cavity views of the plot command files mcnp/plot_* (FIGURES), artemis_<code>_<name>.png:

    neutron_flux_total_<view>         neutron flux
    neutron_flux_<group>_<view>       neutron flux in each energy group (NEUTRON_GROUPS)
    neutron_flux_si_eq_<view>         1-MeV(Si) equivalent neutron flux
    photon_flux_<view>                photon flux
    kerma_neutron_<view>              neutron kerma in silicon
    kerma_photon_<view>               photon kerma in silicon
    kerma_total_<view>                total (neutron + photon) kerma in silicon
    neutron_hardness_<view>           1-MeV(Si) equivalent over total neutron flux
    neutron_fraction_high_<view>      neutron flux above 10 MeV over total neutron flux

Tables of the peak neutron flux in each cassette region by energy group, and of the responses in
silicon in the cassette regions and at their peak, are printed. Reading MCNP results needs
mcnptools (which comes with MCNP), and reading OpenMC results needs the openmc Python package.

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
OpenMC's by the strength of the source. Fluxes are per cm^2 per second, and kerma is in rad(Si)
per second. The responses in silicon (the 1-MeV(Si) equivalent flux and the kerma) come from the
same data in both codes (see the MCNP deck and data/silicon_*.txt).
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

# The quantities tallied: neutron flux (in energy groups), photon flux, 1-MeV(Si) equivalent
# neutron flux, and neutron and photon kerma in silicon
QUANTITIES = ["neutron", "photon", "si", "nkerma", "pkerma"]
# Cell tallies of the responses in silicon in the cassette regions (cells 106 and 206, -x then +x)
CELL_TALLIES = {"si": 34, "nkerma": 44, "pkerma": 54}
# Mesh tallies: each quantity in each plane and view
PLANES = ["xy", "xz"]
VIEWS = ["pit", "cavity"]
MESH_TALLIES = {
    (quantity, plane, view): 64 + 10 * ((k * len(VIEWS) + m) * len(PLANES) + n)
    for k, quantity in enumerate(QUANTITIES)
    for m, view in enumerate(VIEWS)
    for n, plane in enumerate(PLANES)
}
# Peak values in each cassette region: one-bin mesh tallies over a plate as close to the target
# chamber as a part in the cassette can get (from x = +-7.8 cm), at the height of the peak flux,
# -x then +x: neutrons in the four groups, the neutron and photon spectra (same bins as the
# cassette spectra), and the responses in silicon
PEAK_TALLIES = {
    "neutron": {"-x": 264, "+x": 274},
    "si": {"-x": 324, "+x": 334},
    "nkerma": {"-x": 344, "+x": 354},
    "pkerma": {"-x": 364, "+x": 374},
}
PEAK_SPECTRUM_TALLIES = {"neutron": [284, 294], "photon": [304, 314]}
SIDES = ["-x", "+x"]

# Energy groups of the neutron mesh tallies (emesh in the MCNP deck, an energy filter in the
# OpenMC model): name (in the file names) and title of each, in order of energy
NEUTRON_GROUPS = [
    ("thermal", "below 0.5 eV"),
    ("epithermal", "from 0.5 eV to 0.1 MeV"),
    ("fast", "from 0.1 to 10 MeV"),
    ("high", "above 10 MeV"),
]
GROUP = {name: g for g, (name, _) in enumerate(NEUTRON_GROUPS)}

# Map figures: (quantity, energy group (None for all energies), view, color scale (type, vmin,
# vmax; None for the maximum), axis tick spacing [cm]). Quantities: those tallied (QUANTITIES) and
# those derived from them (DERIVED). Log scales: the flux falls by orders of magnitude from the
# source in the chamber to the cassettes (cavity views) and through the water (pit views)
# fmt: off
FIGURES = [
    ('neutron'      , None        , 'pit'   , ('log', 1e6, None) , 100),
    ('neutron'      , None        , 'cavity', ('log', 1e8, None) , 20),
    ('neutron'      , 'thermal'   , 'pit'   , ('log', 1e6, None) , 100),
    ('neutron'      , 'thermal'   , 'cavity', ('log', 1e7, None) , 20),
    ('neutron'      , 'epithermal', 'pit'   , ('log', 1e6, None) , 100),
    ('neutron'      , 'epithermal', 'cavity', ('log', 1e7, None) , 20),
    ('neutron'      , 'fast'      , 'pit'   , ('log', 1e6, None) , 100),
    ('neutron'      , 'fast'      , 'cavity', ('log', 1e7, None) , 20),
    ('neutron'      , 'high'      , 'pit'   , ('log', 1e6, None) , 100),
    ('neutron'      , 'high'      , 'cavity', ('log', 1e7, None) , 20),
    ('si'           , None        , 'pit'   , ('log', 1e6, None) , 100),
    ('si'           , None        , 'cavity', ('log', 1e8, None) , 20),
    ('photon'       , None        , 'pit'   , ('log', 1e6, None) , 100),
    ('photon'       , None        , 'cavity', ('log', 1e8, None) , 20),
    ('nkerma'       , None        , 'pit'   , ('log', 1e-4, None), 100),
    ('nkerma'       , None        , 'cavity', ('log', 1e-1, None), 20),
    ('pkerma'       , None        , 'pit'   , ('log', 1e-4, None), 100),
    ('pkerma'       , None        , 'cavity', ('log', 1e-1, None), 20),
    ('kerma'        , None        , 'pit'   , ('log', 1e-4, None), 100),
    ('kerma'        , None        , 'cavity', ('log', 1e-1, None), 20),
    ('hardness'     , None        , 'pit'   , ('linear', 0, 2)   , 100),
    ('hardness'     , None        , 'cavity', ('linear', 0, 2)   , 20),
    ('high_fraction', None        , 'pit'   , ('linear', 0, 1)   , 100),
    ('high_fraction', None        , 'cavity', ('linear', 0, 1)   , 20),
]
# fmt: on
# Derived quantities: the total kerma (neutron + photon), and the hardness factor and the fraction
# above 10 MeV of the neutron flux. The ratios are left blank where the neutron flux is below the
# floor of the view's flux map, where they are mostly noise
DERIVED = ["kerma", "hardness", "high_fraction"]
RATIO_FLOOR = {"pit": 1e6, "cavity": 1e8}
# Each quantity's title, colorbar label, and file name (with the energy group and view added)
LABELS = {
    "neutron": ("Neutron flux", "n/cm$^2$-s", "neutron_flux"),
    "photon": ("Photon flux", r"$\gamma$/cm$^2$-s", "photon_flux"),
    "si": ("1-MeV(Si) equivalent neutron flux", "n/cm$^2$-s (1-MeV(Si) equivalent)",
           "neutron_flux_si_eq"),
    "nkerma": ("Neutron kerma in silicon", "rad(Si)/s", "kerma_neutron"),
    "pkerma": ("Photon kerma in silicon", "rad(Si)/s", "kerma_photon"),
    "kerma": ("Total kerma in silicon, neutrons and photons", "rad(Si)/s", "kerma_total"),
    "hardness": ("Hardness factor of the neutron flux",
                 "1-MeV(Si) equivalent / total", "neutron_hardness"),
    "high_fraction": ("Fraction of the neutron flux above 10 MeV", "fraction above 10 MeV",
                      "neutron_fraction_high"),
}  # fmt: skip

# Cassette region of the +x cavity (the -x one is its mirror image), as in the models [cm]: a box,
# outside the cavity wall around the target chamber. It is air, like the cavity around it, so the
# flux maps draw its open sides as dashed lines (its other sides are walls)
CASSETTE_X = (1.27, 51.181)
CASSETTE_Y = 20.32  # half-width
CASSETTE_Z = (16.51, 94.615)
CASSETTE_R = 7.06501  # cavity wall around the target chamber: its radius, about x = CASSETTE_XC
CASSETTE_XC = 0.47752
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


def contour_label(x):
    """A contour label: written out from 0.001 up to (not including) 10000, such as 0.3, 1, or 2000,
    and otherwise as a power of ten in the colorbar's style, such as 2×10⁸ (10⁹ for 1×10⁹), in
    mathtext."""
    mantissa, exponent = f"{x:.2e}".split("e")
    mantissa, exponent = f"{float(mantissa):g}", int(exponent)
    if -3 <= exponent < 4:
        return f"{float(mantissa) * 10.0**exponent:.10g}"
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
        """Return (energy edges [MeV], flux in each bin [/cm^2-s] in each cassette region)."""
        t = self.mctal.GetTally(SPECTRUM_TALLIES[particle])
        upper = np.array(t.GetEBins())  # upper edges; the first is the top of a bin from 0
        flux = np.array(
            [
                [t.GetValue(c, 0, 0, 0, 0, 0, e, 0) for e in range(len(upper))]
                for c in range(len(CASSETTE_CELLS))
            ]
        )
        return upper, flux[:, 1:]

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

    def cell(self, quantity):
        """Return (value, relative error) of a cell tally in each cassette region (-x, +x)."""
        t = self.mctal.GetTally(CELL_TALLIES[quantity])
        return [(t.GetValue(c, 0, 0, 0, 0, 0, 0, 0), t.GetError(c, 0, 0, 0, 0, 0, 0, 0))
                for c in range(len(CASSETTE_CELLS))]  # fmt: skip

    def peak(self, quantity, side):
        """Return (value, relative error) of a peak tally in each energy bin and in total (last; the
        only one for a tally with one bin)."""
        t = self.meshtal.GetTally(PEAK_TALLIES[quantity][side])
        n = len(t.GetEBounds()) - 1
        bins = range(n + 1) if n > 1 else range(1)  # the bins, then the total
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
        """Return (energy edges [MeV], flux in each bin [/cm^2-s] in each cassette region)."""
        import openmc

        t = self.statepoint.get_tally(id=SPECTRUM_TALLIES[particle])
        edges = t.find_filter(openmc.EnergyFilter).values * 1e-6
        track = t.mean.reshape(len(CASSETTE_CELLS), -1)
        return edges, track / np.array(self.volumes)[:, None]

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

    def cell(self, quantity):
        """Return (value, relative error) of a cell tally in each cassette region (-x, +x)."""
        t = self.statepoint.get_tally(id=CELL_TALLIES[quantity])
        mean, sd = t.mean.ravel(), t.std_dev.ravel()
        return [(mean[c] / self.volumes[c], sd[c] / mean[c]) for c in range(len(CASSETTE_CELLS))]

    def peak(self, quantity, side):
        """Return (value, relative error) of a peak tally in each energy bin and in total (last; the
        only one for a tally with one bin). The total's error treats the bins as independent."""
        import openmc

        t = self.statepoint.get_tally(id=PEAK_TALLIES[quantity][side])
        mesh = t.find_filter(openmc.MeshFilter).mesh
        volume = np.prod(np.array(mesh.upper_right) - np.array(mesh.lower_left))
        mean, sd = t.mean.ravel() / volume, t.std_dev.ravel() / volume
        if len(mean) > 1:
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
        inner = CASSETTE_XC + np.sqrt(CASSETTE_R**2 - at**2) if abs(at) < CASSETTE_R else x0
        segments = [([inner, x1], [z0, z0]), ([inner, x1], [z1, z1]), ([x1, x1], [z0, z1])]
    return [([sign * x for x in xs], ys) for xs, ys in segments for sign in [1, -1]]


def plot_flux(results):
    """One figure per code and entry of FIGURES, with a panel for each plane."""
    tallies = {}  # each mesh tally, read once

    def panel(r, quantity, group, plane, view):
        """(extent, values) of a quantity, tallied or derived, in a plane."""
        if quantity == "kerma":
            extent, neutron = panel(r, "nkerma", None, plane, view)
            return extent, neutron + panel(r, "pkerma", None, plane, view)[1]
        if quantity in ["hardness", "high_fraction"]:
            top = ("si", None) if quantity == "hardness" else ("neutron", "high")
            extent, values = panel(r, *top, plane, view)
            total = panel(r, "neutron", None, plane, view)[1]
            keep = total >= RATIO_FLOOR[view]
            return extent, np.where(keep, values / np.where(keep, total, 1), np.nan)
        key = r.name, quantity, plane, view
        if key not in tallies:
            tallies[key] = r.mesh(quantity, plane, view)
        return in_plane(*tallies[key], plane, group)

    for quantity, group, view, (scale, vmin, vmax), tick in FIGURES:
        maps = {r.name: [panel(r, quantity, group, plane, view) for plane in PLANES]
                for r in results}  # fmt: skip
        # Where each plane is, across it: the middle of the mesh's single bin across the plane
        across = {}
        for plane in PLANES:
            edges = tallies[results[0].name, "neutron", plane, view][0]
            across[plane] = np.mean(edges[3 - "xyz".index(plane[0]) - "xyz".index(plane[1])][:2])
        title, units, stem = LABELS[quantity]
        if group is not None:
            title += " " + dict(NEUTRON_GROUPS)[group]
        if view == "cavity":
            title += " (cavities)"
        # The neutron flux in total is named like its energy groups
        if group is not None:
            suffix = f"{group}_{view}"
        else:
            suffix = f"total_{view}" if quantity == "neutron" else view
        top = max(np.nanmax(flux) for panels in maps.values() for _, flux in panels)
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
                ax.clabel(cs, fontsize=6, fmt=contour_label if scale == "log" else "%g")
            if levels:
                cbar.add_lines(cs)
            fig.suptitle(f"Artemis: {title}, {name}")
            fig.get_layout_engine().set(h_pad=0.2, wspace=0)
            save(fig, f"artemis_{name.lower()}_{stem}_{suffix}.png", indexed=True)


def plot_spectra(results):
    """The spectra in the cassette regions and at their peak flux (on the plates), per unit
    lethargy, with every code on the same axes, in SHINE colors: blues for the plates, grays for the
    cassette regions, the darker shade for MCNP."""
    for particle in ["neutron", "photon"]:
        fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
        for peak in [True, False]:  # the cassette means (smoother) on top
            for r in results:
                if peak:
                    edges, flux = r.peak_spectrum(particle)
                else:
                    edges, flux = r.spectrum(particle)
                    flux = flux.mean(axis=0)
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
        ax.set_xlabel(f"{particle.capitalize()} energy [MeV]")
        ax.set_ylabel(
            f"{particle.capitalize()} flux per unit lethargy "
            f"[{'n' if particle == 'neutron' else '$\\gamma$'}/cm$^2$-s]"
        )
        ax.set_title(f"Artemis: {particle} spectrum in the cassette regions (mean of -x and +x)")
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
        for side in SIDES:
            flux, error = r.peak("neutron", side)
            cells = "".join(f"{f'{v:.3e} ({e:.3f})':>24}" for v, e in zip(flux, error))
            print(f"{f'{r.name} {side}':14}" + cells)


def silicon_values(r):
    """{(row, where, side): (value, relative error or None)}: the responses in silicon in the
    cassette regions and at their peak, and the quantities derived from them."""
    out = {}
    edges, spectrum = r.spectrum("neutron")
    high = edges[:-1] >= 10.0 - 1e-9  # bins above 10 MeV
    for c, side in enumerate(SIDES):
        cassette = {q: r.cell(q)[c] for q in ["si", "nkerma", "pkerma"]}
        cassette["total"] = (spectrum[c].sum(), None)
        cassette["high"] = (spectrum[c][high].sum(), None)
        groups, errors = r.peak("neutron", side)
        peak = {q: tuple(x[-1] for x in r.peak(q, side)) for q in ["si", "nkerma", "pkerma"]}
        peak["total"] = (groups[-1], errors[-1])
        peak["high"] = (groups[GROUP["high"]], errors[GROUP["high"]])
        for where, v in [("cassette", cassette), ("peak", peak)]:
            out["si", where, side] = v["si"]
            out["nkerma", where, side] = v["nkerma"]
            out["pkerma", where, side] = v["pkerma"]
            out["kerma", where, side] = (v["nkerma"][0] + v["pkerma"][0], None)
            out["hardness", where, side] = (v["si"][0] / v["total"][0], None)
            out["high_fraction", where, side] = (v["high"][0] / v["total"][0], None)
    return out


def print_silicon(results):
    """A table of the responses in silicon in the cassette regions and at their peak, and the
    quantities derived from them."""
    rows = [
        ("si", "1-MeV(Si) equivalent flux [n/cm^2-s]"),
        ("nkerma", "Neutron kerma [rad(Si)/s]"),
        ("pkerma", "Photon kerma [rad(Si)/s]"),
        ("kerma", "Total kerma [rad(Si)/s]"),
        ("hardness", "1-MeV(Si) eq. / total neutron flux"),
        ("high_fraction", "Neutron flux > 10 MeV / total"),
    ]
    columns = [(w, s) for w in ["cassette", "peak"] for s in SIDES]
    print("\nIn silicon, in the cassette regions (mean) and at their peak (relative error)")
    print(f"{'':46}" + "".join(f"{f'{w} {s}':>20}" for w, s in columns))
    values = {r.name: silicon_values(r) for r in results}
    for row, label in rows:
        for r in results:
            cells = []
            for w, s in columns:
                v, e = values[r.name][row, w, s]
                cells.append(f"{f'{v:.3e}' + (f' ({e:.3f})' if e is not None else ''):>20}")
            print(f"{f'{label}, {r.name}':46}" + "".join(cells))


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
    print_silicon(results)


if __name__ == "__main__":
    main()
