#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies, LLC. MIT License: see LICENSE.
# ******************************************************************************

"""Plot the results of the models of each system: the neutron and photon spectra at its peak, and
the neutron and photon flux maps of the mesh tallies, from MCNP, OpenMC, or both.

    python plot_results.py                     every system, both codes
    python plot_results.py --system apollo     one system (SYSTEMS: artemis, apollo)
    python plot_results.py --mcnp              MCNP results (mcnp/<system>/<system>.imsht)
    python plot_results.py --openmc            OpenMC results (openmc/<system>/statepoint.*.h5)

A code whose results are missing is skipped. The plots go in images/<system>/, with the geometry
plots: <system>_spectrum_<particle>.png, the spectra at the peak, with every code given on the same
axes, and maps for the whole and zoomed views of the plot command files mcnp/<system>/plot_* (in
the xy and xz planes; the yz maps are tallied but not plotted), <system>_<code>_<name>.png:

    neutron_flux_total_<view>         neutron flux
    neutron_flux_<group>_<view>       neutron flux in each energy group (NEUTRON_GROUPS)
    neutron_flux_si_eq_<view>         1-MeV(Si) equivalent neutron flux
    photon_flux_<view>                photon flux
    kerma_neutron_<view>              neutron kerma in silicon
    kerma_photon_<view>               photon kerma in silicon
    kerma_total_<view>                total (neutron + photon) kerma in silicon
    neutron_hardness_<view>           1-MeV(Si) equivalent over total neutron flux
    neutron_fraction_high_<view>      neutron flux above 10 MeV over total neutron flux

With --docs, each system's documentation figures (doc_figures) are also copied to docs/figures/.
Tables of the peak neutron flux by energy group, and of the responses in silicon at the peak, are
printed. Reading MCNP results needs mcnptools (which comes with MCNP), and reading OpenMC results
needs the openmc Python package.

The peak is as close to the target chamber as a part can get, at the height of the peak flux: in
Artemis, a 5 x 5 x 0.1 cm plate in each cassette region (-x and +x; one-bin mesh tallies); in
Apollo, a cylindrical shell against the chamber wall, 0.1 cm thick and 5 cm tall (a one-bin
cylindrical mesh tally).

Each flux figure has a panel for each plane of a view (xy and xz), all the same height, with a
shared color scale, the same for every code given: a vmax of None is the largest value in the
panels. All the flux maps of a view (neutrons in total and in each group, 1-MeV(Si) equivalent, and
photons) share one scale, from the view's floor (vmin) to the largest value of them all, and so do
all its dose (kerma) maps (SCALE_GROUP). Anything below vmin is shown in the vmin color instead of
being left blank. There is a contour line at each colorbar tick, labeled with its value and drawn on
the colorbar too. Material boundaries are drawn on top from the OpenMC geometry plot of the same
window (images/<system>/<system>_openmc_geom_<plane>_<view>.png, made by plot_geometry.py), which
shows only material boundaries (the MCNP plots also show the boundaries between cells of one
material). Artemis's cassette regions, which are air like the cavities around them, have their open
sides dashed. Each mesh bin is drawn as a solid block. The figures are about as wide as a GitHub
page (at 150 dpi) and are saved as 256-color indexed PNGs, which look the same as full-color ones
but are much smaller.

All results are absolute: MCNP's are scaled by the total neutron yield with its factor keywords, and
OpenMC's by the strength of the source. Fluxes are per cm^2 per second, and kerma is in rad(Si) per
second. The responses in silicon (the 1-MeV(Si) equivalent flux and the kerma) come from the same
data in both codes (see the MCNP decks and data/silicon_*.txt).
"""

import argparse
import glob
import io
import os
import re
import shutil
import time
from dataclasses import dataclass, field

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS_DIR = os.path.join(HERE, "docs", "figures")

# Tally numbers, the same in every model (see the header of the MCNP decks). The tens digit q of the
# peak tallies is the quantity; 2q4 and 3q4 are the peak values (3q4 only in Artemis, for its +x
# cassette region; Artemis also has cell tallies in its cassette regions, 1q4). The maps take a
# hundred per quantity (4xx to 8xx), with the window in the tens digit: 1 to 3 for the xy, xz, and
# yz planes of the whole view, 4 to 6 for those of the zoomed view.
Q = {"neutron": 1, "nspectrum": 2, "pspectrum": 3, "si": 4, "nkerma": 5, "pkerma": 6}
MAP_BASE = {"neutron": 400, "photon": 500, "si": 600, "nkerma": 700, "pkerma": 800}
WINDOW = {"xy": 1, "xz": 2, "yz": 3}
SPECTRUM_Q = {"neutron": "nspectrum", "photon": "pspectrum"}
# The quantities mapped: neutron flux (in energy groups), photon flux, 1-MeV(Si) equivalent
# neutron flux, and neutron and photon kerma in silicon
QUANTITIES = ["neutron", "photon", "si", "nkerma", "pkerma"]
# The planes plotted (the yz maps are tallied, but not plotted)
PLANES = ["xy", "xz"]

# Energy groups of the neutron mesh tallies (emesh in the MCNP decks, an energy filter in the
# OpenMC models): name (in the file names) and title of each, in order of energy
NEUTRON_GROUPS = [
    ("thermal", "below 0.5 eV"),
    ("epithermal", "from 0.5 eV to 0.1 MeV"),
    ("fast", "from 0.1 to 10 MeV"),
    ("high", "above 10 MeV"),
]
GROUP = {name: g for g, (name, _) in enumerate(NEUTRON_GROUPS)}

# Map figures: (quantity, energy group (None for all energies)), for each view. Quantities: those
# tallied (QUANTITIES) and those derived from them (DERIVED). Fluxes and doses are on log scales
# from each view's floor (System.vmin; the flux falls by orders of magnitude from the source in the
# chamber), the ratios on fixed linear scales (RATIO_SCALES)
MAPS = [
    ("neutron", None),
    ("neutron", "thermal"),
    ("neutron", "epithermal"),
    ("neutron", "fast"),
    ("neutron", "high"),
    ("si", None),
    ("photon", None),
    ("nkerma", None),
    ("pkerma", None),
    ("kerma", None),
    ("hardness", None),
    ("high_fraction", None),
]
# Maps that share one color scale in each view (the view's floor and the largest value of them
# all): the fluxes, and the doses (kerma); the ratios keep their own fixed scales
SCALE_GROUP = {
    "neutron": "flux",
    "photon": "flux",
    "si": "flux",
    "nkerma": "dose",
    "pkerma": "dose",
    "kerma": "dose",
}
RATIO_SCALES = {"hardness": (0, 2), "high_fraction": (0, 1)}
# Derived quantities: the total kerma (neutron + photon), and the hardness factor and the fraction
# above 10 MeV of the neutron flux. The ratios are left blank where the neutron flux is below the
# floor of the view's flux map, where they are mostly noise
DERIVED = ["kerma", "hardness", "high_fraction"]
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
OUTLINE_STYLE = {"color": "dimgray", "linestyle": (0, (4, 2)), "linewidth": 1.2}


# ------------------------------------------------------------------------------
# Systems
# ------------------------------------------------------------------------------
@dataclass
class System:
    """What differs between the systems' models and plots."""

    name: str  # as in the folder and file names
    title: str  # in the plot titles
    views: list  # the whole view, then the zoomed one (as in the plot command files)
    view_titles: dict  # each view's name in the map titles
    regions: list  # where the peak tallies are, such as -x and +x ("" for one)
    peak_hundreds: dict  # hundreds digit of the peak tallies of each region
    mcnp_bin: int  # MCNP radial bin of the peak meshes (1: the first is inside the shell)
    legend: str  # spectrum legend of the peak, after the code
    spectrum_title: str  # spectrum plot title, after the particle
    # (plane, position across it) -> dashed outline segments [(xs, ys)] of a region parts go in,
    # or None
    outline: object
    vmin: dict  # color scale floor of each (scale group, view)
    ratio_floor: dict  # neutron flux below which the ratios are left blank, in each view
    contour_steps: dict  # contour levels in each decade, in each view
    tick: dict  # axis tick spacing [cm] in each view
    height: dict  # panel height [in] in each view
    doc_figures: list = field(default_factory=list)  # figures copied to docs/figures/ by --docs

    @property
    def mcnp_dir(self):
        return os.path.join(HERE, "mcnp", self.name)

    @property
    def openmc_dir(self):
        return os.path.join(HERE, "openmc", self.name)

    @property
    def images_dir(self):
        return os.path.join(HERE, "images", self.name)


def cassette_sides(plane, at):
    """Artemis: the open sides of both cassette regions in a plane at at [cm] across it, as
    (horizontal, vertical) coordinate pairs of line segments. Each cassette region is a box beside
    the cavity wall around the chamber: x = 1.27 to 51.181 cm (mirrored for -x), y = +-20.32 cm,
    z = 16.51 to 94.615 cm; its other sides are walls."""
    (x0, x1), y, (z0, z1) = (1.27, 51.181), 20.32, (16.51, 94.615)
    r, xc = 7.06501, 0.47752  # the cavity wall around the chamber: its radius, about x = xc
    segments = []
    if plane == "xy" and z0 < at < z1:
        # The sides at +-y and the outer end; the inner end is the cavity walls
        segments = [([x0, x1], [y, y]), ([x0, x1], [-y, -y]), ([x1, x1], [-y, y])]
    elif plane == "xz" and abs(at) < y:
        # The bottom, top, and outer end; the inner end is the cavity wall around the chamber
        inner = xc + np.sqrt(r**2 - at**2) if abs(at) < r else x0
        segments = [([inner, x1], [z0, z0]), ([inner, x1], [z1, z1]), ([x1, x1], [z0, z1])]
    return [([sign * x for x in xs], ys) for xs, ys in segments for sign in [1, -1]]


SYSTEMS = {
    "artemis": System(
        name="artemis",
        title="Artemis",
        views=["pit", "cavity"],
        view_titles={"pit": "pit", "cavity": "cavities"},
        regions=["-x", "+x"],
        peak_hundreds={"-x": 200, "+x": 300},
        mcnp_bin=0,
        legend="peak (5 x 5 x 0.1 cm plates, mean of -x and +x)",
        spectrum_title="spectrum at the peak in the cassette regions",
        outline=cassette_sides,
        vmin={("flux", "pit"): 1e6, ("flux", "cavity"): 1e7,
              ("dose", "pit"): 1e-4, ("dose", "cavity"): 1e-1},
        ratio_floor={"pit": 1e6, "cavity": 1e8},
        contour_steps={"pit": [1, 3], "cavity": [1, 2, 5]},
        tick={"pit": 100, "cavity": 20},
        height={"pit": 3.0, "cavity": 3.0},
        doc_figures=[
            "artemis_spectrum_neutron.png",
            "artemis_spectrum_photon.png",
            "artemis_mcnp_neutron_flux_total_cavity.png",
            "artemis_mcnp_neutron_flux_high_cavity.png",
            "artemis_mcnp_neutron_flux_si_eq_cavity.png",
            "artemis_mcnp_neutron_fraction_high_cavity.png",
            "artemis_mcnp_kerma_photon_cavity.png",
        ],
    ),
    "apollo": System(
        name="apollo",
        title="Apollo",
        views=["room", "chamber"],
        view_titles={"room": "room", "chamber": "chamber"},
        regions=[""],
        peak_hundreds={"": 200},
        mcnp_bin=1,
        legend="peak (shell r = 4.9 to 5 cm, 5 cm tall)",
        spectrum_title="spectrum at the peak, against the target chamber",
        outline=None,
        vmin={("flux", "room"): 1e6, ("flux", "chamber"): 1e8,
              ("dose", "room"): 1e-4, ("dose", "chamber"): 1e-1},
        ratio_floor={"room": 1e6, "chamber": 1e8},
        contour_steps={"room": [1, 3], "chamber": [1, 2, 5]},
        tick={"room": 200, "chamber": 20},
        height={"room": 2.2, "chamber": 3.0},
        doc_figures=[
            "apollo_spectrum_neutron.png",
            "apollo_spectrum_photon.png",
            "apollo_mcnp_neutron_flux_total_chamber.png",
            "apollo_mcnp_neutron_flux_high_chamber.png",
            "apollo_mcnp_neutron_flux_si_eq_chamber.png",
            "apollo_mcnp_neutron_fraction_high_chamber.png",
            "apollo_mcnp_kerma_photon_chamber.png",
        ],
    ),
}  # fmt: skip

# SHINE brand colors for the spectrum plots: the blue family and the grays (the product-line
# colors are only for series that are the product lines)
SHINE = {
    "dark blue": "#0093FF",
    "shine blue": "#00CCFF",
    "dark gray": "#53565A",
    "mid gray": "#B1B3B3",
    "light gray": "#CCCCCC",
}
SPECTRUM_COLORS = {"MCNP": SHINE["dark blue"], "OpenMC": SHINE["shine blue"]}

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


def mesh_number(system, quantity, plane, view):
    """Number of a map tally."""
    zoom = system.views.index(view)
    return MAP_BASE[quantity] + 10 * (3 * zoom + WINDOW[plane]) + 4


def peak_number(system, q, region):
    """Number of a peak tally of quantity q (Q) in a region (-x or +x in Artemis)."""
    return system.peak_hundreds[region] + 10 * Q[q] + 4


# ------------------------------------------------------------------------------
# Reading results
# ------------------------------------------------------------------------------
def expand_mcnp(tokens):
    """Numbers of an MCNP list, with its ni and nilog shorthand expanded."""
    out, pending = [], None
    for token in tokens:
        m = re.fullmatch(r"(\d*)i(log)?", token)
        if m:
            pending = (int(m.group(1) or 1), bool(m.group(2)))
            continue
        v = float(token)
        if pending:
            n, log = pending
            a = out[-1]
            for k in range(1, n + 1):
                out.append(a * (v / a) ** (k / (n + 1)) if log else a + (v - a) * k / (n + 1))
            pending = None
        out.append(v)
    return np.array(out)


class MCNPResults:
    """MCNP results: the meshtal (mesh tallies) file of mcnp6 n=<system>.i."""

    name = "MCNP"

    def __init__(self, system):
        import mcnptools

        self.system = system
        base = os.path.join(system.mcnp_dir, system.name)
        if not os.path.exists(base + ".imsht"):
            raise FileNotFoundError(f"no {os.path.relpath(base + '.imsht', HERE)}")
        self.meshtal = mcnptools.Meshtal(base + ".imsht")
        self.deck = open(base + ".i").read()

    def energy_edges(self, n):
        """Upper energy bin edges [MeV] of tally n, from its e card or emesh keyword in the deck
        (the output files round them); the first is the top of a bin from 0."""
        lines = self.deck.split("\n")
        start = next(i for i, x in enumerate(lines) if re.match(rf"\s+(e|fmesh){n}\b", x))
        tokens = None
        for line in lines[start : start + 20]:
            text = line.split("$")[0]
            if tokens is None:
                m = re.match(rf"\s+e{n}\s+(.*)", text) or re.match(r"\s+emesh\s*=\s*(.*)", text)
                if m:
                    tokens = m.group(1).split()
            elif re.match(r"\s{10,}[0-9]", text) and "=" not in text:
                tokens += text.split()
            else:
                break
        return expand_mcnp(tokens)

    def one_bin(self, n):
        """(value, relative error) of a one-bin mesh tally in each energy bin and in total (last;
        the only one for a tally with one bin)."""
        t = self.meshtal.GetTally(n)
        nbin = len(t.GetEBounds()) - 1
        bins = range(nbin + 1) if nbin > 1 else range(1)
        i = self.system.mcnp_bin
        return (np.array([t.GetValue(i, 0, 0, e) for e in bins]),
                np.array([t.GetError(i, 0, 0, e) for e in bins]))  # fmt: skip

    def peak_spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the peak tallies."""
        q = SPECTRUM_Q[particle]
        upper = self.energy_edges(peak_number(self.system, q, self.system.regions[0]))
        flux = [self.one_bin(peak_number(self.system, q, r))[0][: len(upper)]
                for r in self.system.regions]  # fmt: skip
        return upper, np.mean(flux, axis=0)[1:]

    def peak(self, quantity, region):
        """Return (value, relative error) of a peak tally in each energy bin and in total (last; the
        only one for a tally with one bin)."""
        return self.one_bin(peak_number(self.system, quantity, region))

    def mesh(self, quantity, plane, view):
        """Return (edges along x, y, and z [cm], flux[x, y, z, energy group] [/cm^2-s]) of a mesh
        tally."""
        t = self.meshtal.GetTally(mesh_number(self.system, quantity, plane, view))
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
    """OpenMC results: the latest statepoint in openmc/<system>/."""

    name = "OpenMC"

    def __init__(self, system):
        import openmc

        self.system = system
        files = sorted(
            glob.glob(os.path.join(system.openmc_dir, "statepoint.*.h5")),
            key=lambda f: int(f.split(".")[-2]),
        )
        if not files:
            raise FileNotFoundError(f"no statepoint in {system.openmc_dir}")
        self.statepoint = openmc.StatePoint(files[-1])

    def one_bin(self, n):
        """(mean per bin [/cm^2-s], std. dev. per bin, energy edges [MeV] or None) of a one-bin
        mesh tally."""
        import openmc

        t = self.statepoint.get_tally(id=n)
        volume = t.find_filter(openmc.MeshFilter).mesh.volumes.ravel()[0]
        energy = [f for f in t.filters if type(f) is openmc.EnergyFilter]
        edges = energy[0].values * 1e-6 if energy else None
        return t.mean.ravel() / volume, t.std_dev.ravel() / volume, edges

    def peak_spectrum(self, particle):
        """Return (energy edges [MeV], flux in each bin [/cm^2-s]), the mean of the peak tallies."""
        q = SPECTRUM_Q[particle]
        runs = [self.one_bin(peak_number(self.system, q, r)) for r in self.system.regions]
        return runs[0][2], np.mean([mean for mean, _, _ in runs], axis=0)

    def peak(self, quantity, region):
        """Return (value, relative error) of a peak tally in each energy bin and in total (last; the
        only one for a tally with one bin). The total's error treats the bins as independent."""
        mean, sd, _ = self.one_bin(peak_number(self.system, quantity, region))
        if len(mean) > 1:
            mean = np.append(mean, mean.sum())
            sd = np.append(sd, np.sqrt(np.sum(sd**2)))
        return mean, np.divide(sd, mean, out=np.zeros_like(sd), where=mean > 0)

    def mesh(self, quantity, plane, view):
        """Return (edges along x, y, and z [cm], flux[x, y, z, energy group] [/cm^2-s]) of a mesh
        tally."""
        import openmc

        t = self.statepoint.get_tally(id=mesh_number(self.system, quantity, plane, view))
        mesh = t.find_filter(openmc.MeshFilter).mesh
        nx, ny, nz = mesh.dimension
        edges = [
            np.linspace(lo, hi, n + 1)
            for lo, hi, n in zip(mesh.lower_left, mesh.upper_right, mesh.dimension)
        ]
        volume = mesh.volumes.ravel()[0]
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
def geometry_edges(system, plane, view, width=1):
    """Material boundaries from the OpenMC geometry plot of a window, as an RGBA image (black where
    neighboring pixels differ in color), or None if there is no plot. The edges are widened by width
    pixels on each side, so they stay dark when the image is shrunk to the panel size."""
    filename = os.path.join(system.images_dir, f"{system.name}_openmc_geom_{plane}_{view}.png")
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


def plot_flux(system, results):
    """One figure per code, map (MAPS), and view, with a panel for each plane."""
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
            keep = total >= system.ratio_floor[view]
            return extent, np.where(keep, values / np.where(keep, total, 1), np.nan)
        key = r.name, quantity, plane, view
        if key not in tallies:
            tallies[key] = r.mesh(quantity, plane, view)
        return in_plane(*tallies[key], plane, group)

    figures = [(quantity, group, view) for view in system.views for quantity, group in MAPS]
    # Every figure's panels first, for the shared color scales: the view's floor and the largest
    # value of each scale group in each view
    all_maps = [{r.name: [panel(r, quantity, group, plane, view) for plane in PLANES]
                 for r in results} for quantity, group, view in figures]  # fmt: skip
    shared = {}
    for (quantity, _, view), maps in zip(figures, all_maps):
        if quantity in SCALE_GROUP:
            top = max(np.nanmax(flux) for panels in maps.values() for _, flux in panels)
            key = SCALE_GROUP[quantity], view
            shared[key] = max(shared.get(key, 0.0), top)

    for (quantity, group, view), maps in zip(figures, all_maps):
        # Where each plane is, across it: the middle of the mesh's single bin across the plane
        across = {}
        for plane in PLANES:
            edges = tallies[results[0].name, "neutron", plane, view][0]
            across[plane] = np.mean(edges[3 - "xyz".index(plane[0]) - "xyz".index(plane[1])][:2])
        title, units, stem = LABELS[quantity]
        if group is not None:
            title += " " + dict(NEUTRON_GROUPS)[group]
        title += f" ({system.view_titles[view]})"
        # The neutron flux in total is named like its energy groups
        if group is not None:
            suffix = f"{group}_{view}"
        else:
            suffix = f"total_{view}" if quantity == "neutron" else view
        if quantity in SCALE_GROUP:
            scale, vmin, vmax = "log", system.vmin[SCALE_GROUP[quantity], view], None
            top = shared[SCALE_GROUP[quantity], view]
        else:
            scale, (vmin, vmax) = "linear", RATIO_SCALES[quantity]
            top = vmax
        norm = (LogNorm if scale == "log" else Normalize)(vmin, top if vmax is None else vmax)
        for name, panels in maps.items():
            # The panels all have the same height, each with its own aspect ratio
            height = system.height[view]  # [in]
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
                edges = geometry_edges(system, plane, view)
                if edges is not None:
                    ax.imshow(edges, origin="upper", extent=extent)
                for xs, ys in system.outline(plane, across[plane]) if system.outline else []:
                    ax.plot(xs, ys, **OUTLINE_STYLE)
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
                ax.set_title(f"{plane} plane")
                ax.set_xlabel(f"{plane[0]} [cm]")
                ax.set_ylabel(f"{plane[1]} [cm]")
                ax.xaxis.set_major_locator(MultipleLocator(system.tick[view]))
                ax.yaxis.set_major_locator(MultipleLocator(system.tick[view]))
            cbar = fig.colorbar(im, ax=axs, label=units, shrink=0.8, pad=0.01)

            # Contour lines at the view's steps in each decade (not including the ends of the
            # scale), one color per level, the same in every panel, drawn on the colorbar too, and
            # labeled with their values
            if scale == "log":
                decades = range(
                    int(np.floor(np.log10(norm.vmin))), int(np.ceil(np.log10(norm.vmax))) + 1
                )
                levels = [step * 10.0**d for d in decades for step in system.contour_steps[view]]
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
            fig.suptitle(f"{system.title}: {title}, {name}")
            fig.get_layout_engine().set(h_pad=0.2, wspace=0)
            save(system, fig, f"{system.name}_{name.lower()}_{stem}_{suffix}.png", indexed=True)


def plot_spectra(system, results):
    """The spectra at the peak, per unit lethargy, with every code on the same axes, in SHINE blues
    (the darker for MCNP)."""
    for particle in ["neutron", "photon"]:
        fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
        for r in results:
            edges, flux = r.peak_spectrum(particle)
            ax.stairs(
                # Bins with no score are left out (gaps), instead of drops to the axis
                np.where(flux > 0, flux, np.nan) / np.log(edges[1:] / edges[:-1]),
                edges,
                label=f"{r.name}, {system.legend}",
                color=SPECTRUM_COLORS[r.name],
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
        ax.set_title(f"{system.title}: {particle} {system.spectrum_title}")
        ax.grid(which="major", color=SHINE["light gray"], linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(loc="upper left")
        save(system, fig, f"{system.name}_spectrum_{particle}.png")


def row_label(r, region):
    return f"{r.name} {region}".strip()


def print_peak(system, results):
    """A table of the peak neutron flux, by energy group and in total."""
    columns = [title for _, title in NEUTRON_GROUPS] + ["total"]
    print(f"\n{system.title}: peak neutron flux [n/cm^2-s] (relative error)")
    print(f"{'':14}" + "".join(f"{c:>24}" for c in columns))
    for r in results:
        for region in system.regions:
            flux, error = r.peak("neutron", region)
            cells = "".join(f"{f'{v:.3e} ({e:.4f})':>24}" for v, e in zip(flux, error))
            print(f"{row_label(r, region):14}" + cells)


def silicon_values(system, r):
    """{(row, region): (value, relative error or None)}: the neutron flux and the responses in
    silicon at the peak, and the quantities derived from them."""
    out = {}
    for region in system.regions:
        groups, errors = r.peak("neutron", region)
        v = {q: tuple(x[-1] for x in r.peak(q, region)) for q in ["si", "nkerma", "pkerma"]}
        v["total"] = (groups[-1], errors[-1])
        v["high"] = (groups[GROUP["high"]], errors[GROUP["high"]])
        for row in ["total", "high", "si", "nkerma", "pkerma"]:
            out[row, region] = v[row]
        out["kerma", region] = (v["nkerma"][0] + v["pkerma"][0], None)
        out["hardness", region] = (v["si"][0] / v["total"][0], None)
        out["high_fraction", region] = (v["high"][0] / v["total"][0], None)
    return out


def print_silicon(system, results):
    """A table of the neutron flux and the responses in silicon at the peak, and the quantities
    derived from them."""
    rows = [
        ("total", "Neutron flux [n/cm^2-s]"),
        ("high", "Neutron flux > 10 MeV [n/cm^2-s]"),
        ("si", "1-MeV(Si) equivalent flux [n/cm^2-s]"),
        ("nkerma", "Neutron kerma [rad(Si)/s]"),
        ("pkerma", "Photon kerma [rad(Si)/s]"),
        ("kerma", "Total kerma [rad(Si)/s]"),
        ("hardness", "1-MeV(Si) eq. / total neutron flux"),
        ("high_fraction", "Neutron flux > 10 MeV / total"),
    ]
    print(f"\n{system.title}: at the peak (relative error)")
    print(f"{'':46}" + "".join(f"{f'peak {s}'.strip():>24}" for s in system.regions))
    values = {r.name: silicon_values(system, r) for r in results}
    for row, label in rows:
        for r in results:
            cells = []
            for s in system.regions:
                v, e = values[r.name][row, s]
                cells.append(f"{f'{v:.3e}' + (f' ({e:.4f})' if e is not None else ''):>24}")
            print(f"{f'{label}, {r.name}':46}" + "".join(cells))


def save(system, fig, name, indexed=False):
    """Save a figure in images/<system>/, at 150 dpi; indexed: as a 256-color indexed PNG (about a
    quarter the size of a full-color one, and looks the same: the colormap has 256 colors)."""
    from PIL import Image

    os.makedirs(system.images_dir, exist_ok=True)
    filename = os.path.join(system.images_dir, name)
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
    parser.add_argument("--system", choices=list(SYSTEMS), help="one system only (default: all)")
    parser.add_argument("--mcnp", action="store_true", help="MCNP results only")
    parser.add_argument("--openmc", action="store_true", help="OpenMC results only")
    parser.add_argument(
        "--docs", action="store_true", help="also copy the documentation figures to docs/figures/"
    )
    args = parser.parse_args()
    both = not (args.mcnp or args.openmc)

    for system in [SYSTEMS[args.system]] if args.system else SYSTEMS.values():
        results = []
        for wanted, code in [
            (args.mcnp or both, MCNPResults),
            (args.openmc or both, OpenMCResults),
        ]:
            if wanted:
                try:
                    results.append(code(system))
                except (FileNotFoundError, OSError) as error:
                    print(f"{system.title}: no {code.name} results ({error}); skipped")
        if not results:
            continue
        plot_spectra(system, results)
        plot_flux(system, results)
        print_peak(system, results)
        print_silicon(system, results)
        if args.docs:
            os.makedirs(DOCS_DIR, exist_ok=True)
            for name in system.doc_figures:
                shutil.copy2(os.path.join(system.images_dir, name), os.path.join(DOCS_DIR, name))
                print("Copied", os.path.relpath(os.path.join(DOCS_DIR, name), HERE))


if __name__ == "__main__":
    main()
