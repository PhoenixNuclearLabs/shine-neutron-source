#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies. All rights reserved.
# ******************************************************************************

# A very simplified model of the Artemis neutron generator, the same as ../mcnp/artemis.i: the
# same dimensions, cells, surfaces, materials, and tallies, with the same cell, surface, material,
# and tally numbers. It holds the target chamber with its double wall, the upstream section of the
# beamline (the beam tube runs up to the ground), the two irradiation cavities with their cassette
# regions and ducts, the water pit (full of water up to the ground), its liner, the concrete around
# it, the earth, and the air above ground. The model is symmetric about x = 0 and y = 0. The
# backstop is at the origin and the beam travels down (-z); the ground is at z = 289.56 cm.
#
# Differences from the MCNP model:
# - There is no cell for the void outside the model (MCNP cell 999): the outer surfaces of the
#   model are vacuum boundaries instead.
# - Tallies are already absolute rates: OpenMC scales them by the total strength of the sources,
#   which beam-target gives in n/s (their sum is the total neutron yield), as the MCNP model does
#   with its fm cards. But flux tallies are not divided by volume: divide the cassette tallies by
#   the volume of cells 106 and 206 (the volume of each cell MCNP needs a vol= for is set on the
#   cell) and the mesh tallies by the volume of a mesh element.
#
# The neutron source comes from ../source/Artemis/source.h5, made by beam-target
# (https://github.com/lucasjacobson/beam-target), which must be installed to run this script. The
# geometry plots and mesh tallies follow the MCNP plot command files ../mcnp/plot_*. Set
# OPENMC_CROSS_SECTIONS to the cross_sections.xml of an ENDF/B-VIII.1 library, the data the MCNP
# model uses.
#
# Run it from anywhere; OpenMC's files (model.xml, statepoints, ...) go in this folder:
#   python artemis.py               write model.xml and run (a statepoint after every batch)
#   python artemis.py --restart     continue the run from its latest statepoint
#   python artemis.py --plot        make the geometry plots (../images/artemis_openmc_geom_*.png)
#   python artemis.py --inputs      only write model.xml
# Add --threads N to use N threads (default: all).

import argparse
import functools
import glob
import math
import os
import re

import h5py
import matplotlib.colors
import numpy as np
from beamtarget.openmc_source import get_openmc_sources

import openmc

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE_H5 = os.path.join(HERE, "..", "source", "Artemis", "source.h5")
MCNP_DIR = os.path.join(HERE, "..", "mcnp")
IMAGES_DIR = os.path.join(HERE, "..", "images")

# ------------------------------------------------------------------------------
# Dimensions [cm]
# ------------------------------------------------------------------------------
# fmt: off
GROUND = 289.56
Z_ROOF = round(GROUND + 462.28, 5)       # top of the air above ground

# Target chamber and upstream section: radii
R_BEAM_IN, R_BEAM_OUT = 3.9624, 4.5974   # beam tube
R_BORE = 5.4483                          # target chamber inner radius
R_WALL_1, R_WATER, R_WALL_2 = 5.7277, 5.8674, 6.096
R_GAP, R_CAV_WALL = 6.40969, 7.06501     # water around the chamber, cavity wall around it
R_UP_IN, R_UP_OUT = 12.065, 12.7         # upstream tube

# Target chamber and upstream section: heights
Z_END = -3.52552                         # bottom of the target chamber end plate
Z_END_TOP = -1.18872                     # top of the end plate
Z_PLATE = -0.9525                        # bottom of the steel plate below the beam stop
Z_STOP = -0.635                          # bottom of the beam stop
Z_WALLS = 129.54                         # top of the target chamber walls
Z_UP0, Z_UP1, Z_UP2 = 130.175, 208.026, 208.661  # upstream tube: bottom plate, top, top plate

# Cavities and ducts (+x and +y values, mirrored)
XC0, XC1, XCS, XCI, XCO = 0.0, 1.27, 51.181, 66.675, 67.945
WALL = 1.27                              # cavity and duct walls
XVI0, XVI1 = 128.27, 189.23              # vertical duct, inside
XVO0, XVO1 = XVI0 - WALL, XVI1 + WALL
YCS, YDI, YCI, YCO = 20.32, 22.5425, 44.958, 46.228
YDO = YDI + WALL
Z_CAV0, Z_CAV1 = -4.2672, 103.6828       # cavities, outside: bottom and top
Z_AIR0, Z_AIR1 = -2.9972, 102.4128       # cavities, inside: bottom and top
Z_HD0 = 10.9728                          # horizontal ducts, inside: bottom
Z_HD0_OUT = Z_HD0 - WALL
Z_CAS0, Z_CAS1 = 16.51, 94.615           # cassette regions

# Pit, liner, air gap, concrete, earth
XP, YP = 195.2371, 195.2371
XL, YL = 195.58, 195.58
XA, YA = 215.9, 215.9
XW, YW = 256.54, 256.54
EARTH_T = 3 * 30.48                      # earth: 3 ft thick beside and below the concrete
XE, YE = XW + EARTH_T, YW + EARTH_T      # floor and earth
Z_SLAB, Z_GAP, Z_LINER, Z_PIT, Z_FLOOR = -208.28, -177.8, -167.64, -167.2971, 269.24
Z_EARTH = round(Z_SLAB - EARTH_T, 5)
# fmt: on

# ------------------------------------------------------------------------------
# Materials from the PNNL compendium of material compositions (PNNL-15870), as in the MCNP model:
# number -> (name, density [g/cm3], thermal scattering data, atom fractions)
# ------------------------------------------------------------------------------
COMPENDIUM = {}

# Air (dry, near sea level), rho = 0.001205 g/cm^3
# fmt: off
COMPENDIUM[4] = ('Air (dry, near sea level)', 0.001205, None, [
    ('C12', 1.4858312730e-04), ('C13', 1.6070347337e-06), ('N14', 7.8157536121e-01),
    ('N15', 2.8553277076e-03), ('O16', 2.1023586223e-01), ('O17', 8.0084232333e-05),
    ('O18', 4.3203335864e-04), ('Ar36', 1.5718389806e-05), ('Ar38', 2.9521611760e-06),
    ('Ar40', 4.6524705503e-03),
])
# Aluminum, alloy 6061-O, rho = 2.7 g/cm^3
COMPENDIUM[13] = ('Aluminum, alloy 6061-O', 2.7, 'c_Al27', [
    ('Mg24', 8.8169277578e-03), ('Mg25', 1.1162080969e-03), ('Mg26', 1.2289451147e-03),
    ('Al27', 9.7733008842e-01), ('Si28', 5.3450290793e-03), ('Si29', 2.7153162700e-04),
    ('Si30', 1.7920507805e-04), ('Ti46', 4.0960555374e-05), ('Ti47', 3.6938973573e-05),
    ('Ti48', 3.6601359299e-04), ('Ti49', 2.6860194494e-05), ('Ti50', 2.5718263859e-05),
    ('Cr50', 4.4207409019e-05), ('Cr52', 8.5249587900e-04), ('Cr53', 9.6666189433e-05),
    ('Cr54', 2.4062260605e-05), ('Mn55', 4.3258616015e-04), ('Fe54', 1.1607847426e-04),
    ('Fe56', 1.8221838027e-03), ('Fe57', 4.2082170564e-05), ('Fe58', 5.6003643696e-06),
    ('Cu63', 8.1185431728e-04), ('Cu65', 3.6219386389e-04), ('Zn64', 2.9229283659e-04),
    ('Zn66', 1.6940606828e-04), ('Zn67', 2.4840167724e-05), ('Zn68', 1.1520218205e-04),
    ('Zn70', 3.8210984480e-06),
])
# Concrete, Portland, rho = 2.3 g/cm^3
COMPENDIUM[108] = ('Concrete, Portland', 2.3, None, [
    ('H1', 1.6873924210e-01), ('H2', 1.9407244675e-05), ('C12', 1.4010687221e-03),
    ('C13', 1.5153578617e-05), ('O16', 5.6115663647e-01), ('O17', 2.1375895612e-04),
    ('O18', 1.1531733159e-03), ('Na23', 1.1838220257e-02), ('Mg24', 1.1056218375e-03),
    ('Mg25', 1.3996984903e-04), ('Mg26', 1.5410680379e-04), ('Al27', 2.1353800208e-02),
    ('Si28', 1.8824138727e-01), ('Si29', 9.5628086201e-03), ('Si30', 6.3112495738e-03),
    ('K39', 5.2744001423e-03), ('K40', 6.6171712339e-07), ('K41', 3.8064004990e-04),
    ('Ca40', 1.8103182585e-02), ('Ca42', 1.2082358479e-04), ('Ca43', 2.5210485234e-05),
    ('Ca44', 3.8954868294e-04), ('Ca46', 7.4697734026e-07), ('Ca48', 3.4921190657e-05),
    ('Fe54', 2.4924598409e-04), ('Fe56', 3.9126289178e-03), ('Fe57', 9.0359664721e-05),
    ('Fe58', 1.2025212577e-05),
])
# Copper (Cu), rho = 8.96 g/cm^3
COMPENDIUM[112] = ('Copper (Cu)', 8.96, None, [
    ('Cu63', 6.9150000000e-01), ('Cu65', 3.0850000000e-01),
])
# Earth, U.S. Average, rho = 1.52 g/cm^3
COMPENDIUM[116] = ('Earth, U.S. Average', 1.52, None, [
    ('O16', 6.6897625085e-01), ('O17', 2.5483021274e-04), ('O18', 1.3747419371e-03),
    ('Na23', 5.5780533239e-03), ('Mg24', 9.0303941357e-03), ('Mg25', 1.1432325783e-03),
    ('Mg26', 1.2586990687e-03), ('Al27', 5.3073574433e-02), ('Si28', 1.8597889604e-01),
    ('Si29', 9.4478723089e-03), ('Si30', 6.2353940617e-03), ('K39', 7.1372817329e-03),
    ('K40', 8.9543102717e-07), ('K41', 5.1507947855e-04), ('Ca40', 2.5848275793e-02),
    ('Ca42', 1.7251559647e-04), ('Ca43', 3.5996299110e-05), ('Ca44', 5.5620948106e-04),
    ('Ca46', 1.0665570107e-06), ('Ca48', 4.9861540248e-05), ('Ti46', 1.6573118377e-04),
    ('Ti47', 1.4945939482e-04), ('Ti48', 1.4809336809e-03), ('Ti49', 1.0867947930e-04),
    ('Ti50', 1.0405909478e-04), ('Mn55', 2.7182020995e-04), ('Fe54', 1.2303839620e-03),
    ('Fe56', 1.9314396929e-02), ('Fe57', 4.4605365534e-04), ('Fe58', 5.9361553000e-05),
])
# Steel, Stainless 304, rho = 8.03 g/cm^3
COMPENDIUM[331] = ('Steel, Stainless 304', 8.03, 'c_Fe56', [
    ('C12', 3.5965407501e-03), ('C13', 3.8899207546e-05), ('Si28', 1.7922298752e-02),
    ('Si29', 9.1046669109e-04), ('Si30', 6.0088858247e-04), ('P31', 7.9296709964e-04),
    ('S32', 4.8507331818e-04), ('S33', 3.8299293466e-06), ('S34', 2.1702932964e-05),
    ('S36', 5.1065724621e-08), ('Cr50', 8.6658055923e-03), ('Cr52', 1.6711143493e-01),
    ('Cr53', 1.8949095266e-02), ('Cr54', 4.7168308920e-03), ('Mn55', 1.9869795658e-02),
    ('Fe54', 3.9042933613e-02), ('Fe56', 6.1289056129e-01), ('Fe57', 1.4154315881e-02),
    ('Fe58', 1.8836796029e-03), ('Ni58', 6.0141059320e-02), ('Ni60', 2.3166228378e-02),
    ('Ni61', 1.0070199072e-03), ('Ni62', 3.2108201181e-03), ('Ni64', 8.1770122474e-04),
])
# Water, Liquid (H2O), rho = 0.997 g/cm^3
COMPENDIUM[392] = ('Water, Liquid (H2O)', 0.997, 'c_H_in_H2O', [
    ('H1', 6.6659000000e-01), ('H2', 7.6666666667e-05), ('O16', 3.3252333333e-01),
    ('O17', 1.2666666667e-04), ('O18', 6.8333333333e-04),
])
# fmt: on


def check_source(fname_h5):
    """Check that the source was made for this target chamber."""
    with h5py.File(fname_h5, "r") as hf:
        assert hf["inputs"]["wall_radius"].attrs["units"] == "cm"
        radius = float(hf["inputs"]["wall_radius"][0])
        assert np.isclose(radius, R_BORE), f"{fname_h5} has wall_radius {radius} cm"
        assert np.allclose(hf.attrs["backstop_position"], 0.0)


def get_materials():
    """Materials by MCNP material number."""
    mats = {}
    for number, (name, density, sab, nuclides) in COMPENDIUM.items():
        mat = openmc.Material(number, name)
        for nuclide, fraction in nuclides:
            mat.add_nuclide(nuclide, fraction, "ao")
        mat.set_density("g/cm3", density)
        if sab is not None:
            mat.add_s_alpha_beta(sab)
        mats[number] = mat
    return mats


# ------------------------------------------------------------------------------
# Surfaces, numbered as in the MCNP model: cylinders along z 101+, planes normal to x 201+,
# y 301+, and z 401+, each sorted by position
# ------------------------------------------------------------------------------
def get_surfaces():
    x_planes = [XC1, XCS, XCI, XCO, XVO0, XVI0, XVI1, XVO1, XP, XL, XA, XW, XE]
    y_planes = [YCS, YDI, YDO, YCI, YCO, YP, YL, YA, YW, YE]
    groups = [
        (
            "cz",
            101,
            [
                R_BEAM_IN,
                R_BEAM_OUT,
                R_BORE,
                R_WALL_1,
                R_WATER,
                R_WALL_2,
                R_GAP,
                R_CAV_WALL,
                R_UP_IN,
                R_UP_OUT,
            ],
        ),
        ("px", 201, [XC0] + x_planes + [-v for v in x_planes]),
        ("py", 301, y_planes + [-v for v in y_planes]),
        (
            "pz",
            401,
            [
                Z_EARTH,
                Z_SLAB,
                Z_GAP,
                Z_LINER,
                Z_PIT,
                Z_CAV0,
                Z_END,
                Z_AIR0,
                Z_END_TOP,
                Z_PLATE,
                Z_STOP,
                0.0,
                Z_HD0_OUT,
                Z_HD0,
                Z_CAS0,
                Z_CAS1,
                Z_AIR1,
                Z_CAV1,
                Z_WALLS,
                Z_UP0,
                Z_UP1,
                Z_UP2,
                Z_FLOOR,
                GROUND,
                Z_ROOF,
            ],
        ),
    ]
    make = {
        "cz": lambda n, v: openmc.ZCylinder(surface_id=n, r=v),
        "px": lambda n, v: openmc.XPlane(surface_id=n, x0=v),
        "py": lambda n, v: openmc.YPlane(surface_id=n, y0=v),
        "pz": lambda n, v: openmc.ZPlane(surface_id=n, z0=v),
    }
    surfaces = {}
    for kind, start, values in groups:
        for n, v in enumerate(sorted(set(values)), start=start):
            surfaces[kind, v] = make[kind](n, v)
    # Outer surfaces of the model: vacuum boundaries (MCNP: the importance-0 cell 999 outside)
    for key in [("px", -XE), ("px", XE), ("py", -YE), ("py", YE), ("pz", Z_EARTH), ("pz", Z_ROOF)]:
        surfaces[key].boundary_type = "vacuum"
    return surfaces


# Material numbers (None: void)
VOID, CU, SS, H2O, AIR, AL, CONC, EARTH = None, 112, 331, 392, 4, 13, 108, 116

# Target chamber and upstream section: cylindrical shells (inside r1 only, if r0 is None)
# fmt: off
CHAMBER_CELLS = [
    ( 1, VOID, None,       R_BORE,     0.0,       Z_UP0,     'Vacuum in target chamber'),
    ( 2, CU,   None,       R_BORE,     Z_STOP,    0.0,       'Copper beam stop'),
    ( 3, SS,   None,       R_WALL_1,   Z_PLATE,   Z_STOP,    'Steel plate below beam stop'),
    ( 4, SS,   R_BORE,     R_WALL_1,   Z_STOP,    Z_WALLS,   'Target chamber inner wall'),
    ( 5, H2O,  R_WALL_1,   R_WATER,    Z_PLATE,   Z_WALLS,   'Cooling water jacket'),
    ( 6, H2O,  None,       R_WATER,    Z_END_TOP, Z_PLATE,   'Cooling water behind beam stop'),
    ( 7, SS,   R_WATER,    R_WALL_2,   Z_END_TOP, Z_WALLS,   'Target chamber outer wall'),
    ( 8, SS,   None,       R_WALL_2,   Z_END,     Z_END_TOP, 'Target chamber end plate'),
    ( 9, SS,   R_BORE,     R_UP_OUT,   Z_WALLS,   Z_UP0,     'Upstream tube bottom plate'),
    (10, VOID, None,       R_UP_IN,    Z_UP0,     Z_UP1,     'Vacuum in upstream tube'),
    (11, SS,   R_UP_IN,    R_UP_OUT,   Z_UP0,     Z_UP1,     'Upstream tube'),
    (12, SS,   R_BEAM_OUT, R_UP_OUT,   Z_UP1,     Z_UP2,     'Upstream tube top plate'),
    (13, VOID, None,       R_BEAM_IN,  Z_UP1,     GROUND,    'Vacuum in beam tube'),
    (14, SS,   R_BEAM_IN,  R_BEAM_OUT, Z_UP1,     GROUND,    'Beam tube'),
    (15, H2O,  R_WALL_2,   R_GAP,      Z_END,     Z_CAV1,    'Water between chamber and cavities'),
    (16, H2O,  None,       R_GAP,      Z_CAV0,    Z_END,     'Water below chamber, in cavities'),
]
# fmt: on


# ------------------------------------------------------------------------------
# Cells, numbered as in the MCNP model; each is a plain intersection of half-spaces
# ------------------------------------------------------------------------------
def get_geometry(mats):

    S = get_surfaces()
    cells = []

    def box(x0, x1, y0, y1, z0, z1):
        return (
            +S["px", x0] & -S["px", x1] & +S["py", y0] & -S["py", y1] & +S["pz", z0] & -S["pz", z1]
        )

    def tube(r0, r1, z0, z1):
        """Between two cylinders (r0 = None: inside r1), from z0 to z1."""
        region = -S["cz", r1] & +S["pz", z0] & -S["pz", z1]
        return region if r0 is None else +S["cz", r0] & region

    def cell(number, mat, region, name):
        fill = None if mat is None else mats[mat]
        cells.append(openmc.Cell(cell_id=number, name=name, fill=fill, region=region))

    # Target chamber and upstream section
    for number, mat, r0, r1, z0, z1, name in CHAMBER_CELLS:
        cell(number, mat, tube(r0, r1, z0, z1), name)

    # Irradiation cavities, cassette regions, and ducts, with their walls split into slabs: boxes
    # given for the +x side and mirrored for the -x side
    def side(sign, base, tag):
        def b(x0, x1, y0, y1, z0, z1):
            if sign < 0:
                x0, x1 = -x1, -x0
            return box(x0, x1, y0, y1, z0, z1)

        def cl(i, mat, region, name):
            cell(base + i, mat, region, f"{name} ({tag})")

        out_gap, out_wall = +S["cz", R_GAP], +S["cz", R_CAV_WALL]
        inner_x = +S["px", XC1] if sign > 0 else -S["px", -XC1]
        # fmt: off
        cl( 1, AIR, b(XCS, XCI, -YCI, YCI, Z_AIR0, Z_AIR1), 'Air in cavity, beyond cassette')
        cl( 2, AIR, b(XC1, XCS, YCS, YCI, Z_AIR0, Z_AIR1), 'Air in cavity, +y of cassette')
        cl( 3, AIR, b(XC1, XCS, -YCI, -YCS, Z_AIR0, Z_AIR1), 'Air in cavity, -y of cassette')
        cl( 4, AIR, b(XC1, XCS, -YCS, YCS, Z_AIR0, Z_CAS0) & out_wall,
            'Air in cavity, below cassette')
        cl( 5, AIR, b(XC1, XCS, -YCS, YCS, Z_CAS1, Z_AIR1) & out_wall,
            'Air in cavity, above cassette')
        cl( 6, AIR, b(XC1, XCS, -YCS, YCS, Z_CAS0, Z_CAS1) & out_wall, 'Cassette region')
        cl( 7, AL,  b(XC0, XCO, -YCO, YCO, Z_CAV0, Z_AIR0) & out_gap, 'Cavity bottom')
        cl( 8, AL,  b(XC0, XCO, -YCO, YCO, Z_AIR1, Z_CAV1) & out_gap, 'Cavity top')
        cl( 9, AL,  b(XC0, XCO, YCI, YCO, Z_AIR0, Z_AIR1), 'Cavity +y wall')
        cl(10, AL,  b(XC0, XCO, -YCO, -YCI, Z_AIR0, Z_AIR1), 'Cavity -y wall')
        cl(11, AL,  b(XC0, XC1, -YCI, YCI, Z_AIR0, Z_AIR1) & out_gap, 'Cavity inner wall')
        cl(12, AL,  inner_x & tube(R_GAP, R_CAV_WALL, Z_AIR0, Z_AIR1), 'Cavity wall around chamber')
        cl(13, AL,  b(XCI, XCO, YDI, YCI, Z_AIR0, Z_AIR1), 'Cavity outer wall, +y of duct')
        cl(14, AL,  b(XCI, XCO, -YCI, -YDI, Z_AIR0, Z_AIR1), 'Cavity outer wall, -y of duct')
        cl(15, AL,  b(XCI, XCO, -YDI, YDI, Z_AIR0, Z_HD0), 'Cavity outer wall, below duct')
        cl(16, AIR, b(XCI, XVI0, -YDI, YDI, Z_HD0, Z_AIR1), 'Air in horizontal duct')
        cl(17, AL,  b(XCO, XVO0, YDI, YDO, Z_HD0_OUT, Z_CAV1), 'Horizontal duct +y wall')
        cl(18, AL,  b(XCO, XVO0, -YDO, -YDI, Z_HD0_OUT, Z_CAV1), 'Horizontal duct -y wall')
        cl(19, AL,  b(XCO, XVO0, -YDI, YDI, Z_HD0_OUT, Z_HD0), 'Horizontal duct bottom')
        cl(20, AL,  b(XCO, XVO0, -YDI, YDI, Z_AIR1, Z_CAV1), 'Horizontal duct top')
        cl(21, AIR, b(XVI0, XVI1, -YDI, YDI, Z_AIR0, GROUND), 'Air in vertical duct')
        cl(22, AL,  b(XVO0, XVI0, -YDI, YDI, Z_AIR0, Z_HD0), 'Vertical duct inner wall, below')
        cl(23, AL,  b(XVO0, XVI0, -YDI, YDI, Z_AIR1, GROUND), 'Vertical duct inner wall, above')
        cl(24, AL,  b(XVI1, XVO1, -YDI, YDI, Z_AIR0, GROUND), 'Vertical duct outer wall')
        cl(25, AL,  b(XVO0, XVO1, YDI, YDO, Z_CAV0, GROUND), 'Vertical duct +y wall')
        cl(26, AL,  b(XVO0, XVO1, -YDO, -YDI, Z_CAV0, GROUND), 'Vertical duct -y wall')
        cl(27, AL,  b(XVO0, XVO1, -YDI, YDI, Z_CAV0, Z_AIR0), 'Vertical duct bottom')
        # fmt: on

    side(-1, 100, "-x")
    side(+1, 200, "+x")

    # Water in the pit, up to the ground
    water = [
        (box(-XP, XP, -YP, YP, Z_PIT, Z_CAV0), "Water below everything"),
        (box(-XP, XP, -YP, -YCO, Z_CAV0, GROUND), "Water, -y side"),
        (box(-XP, XP, YCO, YP, Z_CAV0, GROUND), "Water, +y side"),
        (box(-XP, -XVO1, -YCO, YCO, Z_CAV0, GROUND), "Water, -x end"),
        (box(XVO1, XP, -YCO, YCO, Z_CAV0, GROUND), "Water, +x end"),
        (box(-XVO1, -XVO0, -YCO, -YDO, Z_CAV0, GROUND), "Water, -y of vertical duct (-x)"),
        (box(-XVO1, -XVO0, YDO, YCO, Z_CAV0, GROUND), "Water, +y of vertical duct (-x)"),
        (box(XVO0, XVO1, -YCO, -YDO, Z_CAV0, GROUND), "Water, -y of vertical duct (+x)"),
        (box(XVO0, XVO1, YDO, YCO, Z_CAV0, GROUND), "Water, +y of vertical duct (+x)"),
        (box(-XVO0, -XCO, -YCO, -YDO, Z_CAV0, Z_CAV1), "Water, -y of horizontal duct (-x)"),
        (box(-XVO0, -XCO, YDO, YCO, Z_CAV0, Z_CAV1), "Water, +y of horizontal duct (-x)"),
        (box(-XVO0, -XCO, -YDO, YDO, Z_CAV0, Z_HD0_OUT), "Water below horizontal duct (-x)"),
        (box(XCO, XVO0, -YCO, -YDO, Z_CAV0, Z_CAV1), "Water, -y of horizontal duct (+x)"),
        (box(XCO, XVO0, YDO, YCO, Z_CAV0, Z_CAV1), "Water, +y of horizontal duct (+x)"),
        (box(XCO, XVO0, -YDO, YDO, Z_CAV0, Z_HD0_OUT), "Water below horizontal duct (+x)"),
        (
            box(-XVO0, XVO0, -YCO, YCO, Z_CAV1, Z_WALLS) & +S["cz", R_WALL_2],
            "Water around target chamber",
        ),
        (
            box(-XVO0, XVO0, -YCO, YCO, Z_WALLS, Z_UP2) & +S["cz", R_UP_OUT],
            "Water around upstream tube",
        ),
        (
            box(-XVO0, XVO0, -YCO, YCO, Z_UP2, GROUND) & +S["cz", R_BEAM_OUT],
            "Water around beam tube",
        ),
    ]
    for n, (region, name) in enumerate(water, start=301):
        cell(n, H2O, region, name)

    # Pit liner and air gap
    cell(401, SS, box(-XL, XL, -YL, YL, Z_LINER, Z_PIT), "Pit liner, bottom")
    cell(402, SS, box(-XL, XL, -YL, -YP, Z_PIT, GROUND), "Pit liner, -y side")
    cell(403, SS, box(-XL, XL, YP, YL, Z_PIT, GROUND), "Pit liner, +y side")
    cell(404, SS, box(-XL, -XP, -YP, YP, Z_PIT, GROUND), "Pit liner, -x side")
    cell(405, SS, box(XP, XL, -YP, YP, Z_PIT, GROUND), "Pit liner, +x side")
    cell(411, AIR, box(-XA, XA, -YA, YA, Z_GAP, Z_LINER), "Air gap, bottom")
    cell(412, AIR, box(-XA, XA, -YA, -YL, Z_LINER, GROUND), "Air gap, -y side")
    cell(413, AIR, box(-XA, XA, YL, YA, Z_LINER, GROUND), "Air gap, +y side")
    cell(414, AIR, box(-XA, -XL, -YL, YL, Z_LINER, GROUND), "Air gap, -x side")
    cell(415, AIR, box(XL, XA, -YL, YL, Z_LINER, GROUND), "Air gap, +x side")

    # Concrete
    cell(501, CONC, box(-XW, XW, -YW, YW, Z_SLAB, Z_GAP), "Slab below pit")
    cell(502, CONC, box(-XW, XW, -YW, -YA, Z_GAP, Z_FLOOR), "Pit wall, -y side")
    cell(503, CONC, box(-XW, XW, YA, YW, Z_GAP, Z_FLOOR), "Pit wall, +y side")
    cell(504, CONC, box(-XW, -XA, -YA, YA, Z_GAP, Z_FLOOR), "Pit wall, -x side")
    cell(505, CONC, box(XA, XW, -YA, YA, Z_GAP, Z_FLOOR), "Pit wall, +x side")
    cell(506, CONC, box(-XE, XE, -YE, -YA, Z_FLOOR, GROUND), "Floor, -y side")
    cell(507, CONC, box(-XE, XE, YA, YE, Z_FLOOR, GROUND), "Floor, +y side")
    cell(508, CONC, box(-XE, -XA, -YA, YA, Z_FLOOR, GROUND), "Floor, -x side")
    cell(509, CONC, box(XA, XE, -YA, YA, Z_FLOOR, GROUND), "Floor, +x side")

    # Earth
    cell(601, EARTH, box(-XE, XE, -YE, YE, Z_EARTH, Z_SLAB), "Earth below slab")
    cell(602, EARTH, box(-XE, XE, -YE, -YW, Z_SLAB, Z_FLOOR), "Earth beside pit walls, -y side")
    cell(603, EARTH, box(-XE, XE, YW, YE, Z_SLAB, Z_FLOOR), "Earth beside pit walls, +y side")
    cell(604, EARTH, box(-XE, -XW, -YW, YW, Z_SLAB, Z_FLOOR), "Earth beside pit walls, -x side")
    cell(605, EARTH, box(XW, XE, -YW, YW, Z_SLAB, Z_FLOOR), "Earth beside pit walls, +x side")

    # Air above ground
    cell(701, AIR, box(-XE, XE, -YE, YE, GROUND, Z_ROOF), "Air above ground")

    # The volume of each cell MCNP needs a vol= for, as in the MCNP model
    by_id = {surface.id: surface for surface in S.values()}
    for c in cells:
        volume = prism_volume(c.region, by_id)
        if volume is not None:
            c.volume = volume

    return openmc.Geometry(openmc.Universe(cells=cells))


# ------------------------------------------------------------------------------
# Volumes of the cells bounded by both planes normal to x or y and cylinders along z, which MCNP
# can't calculate (they are neither polyhedra nor symmetric about one axis). Each is a prism along
# z, so its volume is its cross-section (a rectangle, inside one cylinder and/or outside another)
# times its height. The cross-section is the integral over x of the chord length in y, in closed
# form piece by piece.
# ------------------------------------------------------------------------------
def disk_rect_area(r, x0, x1, y0, y1):
    """Area of the rectangle [x0, x1] x [y0, y1] inside the disk of radius r about the origin."""
    x0, x1 = max(x0, -r), min(x1, r)
    if x1 <= x0:
        return 0.0
    cuts = {x0, x1}
    for yy in (y0, y1):
        if abs(yy) < r:
            w = math.sqrt(r * r - yy * yy)
            cuts |= {-w, w}
    xs = sorted(c for c in cuts if x0 <= c <= x1)

    def chord(t):  # integral of sqrt(r^2 - x^2) from 0 to t
        t = max(-r, min(r, t))
        return 0.5 * (t * math.sqrt(r * r - t * t) + r * r * math.asin(t / r))

    area = 0.0
    for a, b in zip(xs, xs[1:]):
        m = 0.5 * (a + b)
        h = math.sqrt(r * r - m * m)
        top = chord(b) - chord(a) if h < y1 else y1 * (b - a)
        bottom = -(chord(b) - chord(a)) if -h > y0 else y0 * (b - a)
        if min(h, y1) > max(-h, y0):
            area += top - bottom
    return area


def prism_volume(region, surfaces):
    """Volume of a cell bounded by planes normal to x or y and cylinders along z (and planes
    normal to z), or None for any other cell."""
    lo = {"x": -math.inf, "y": -math.inf, "z": -math.inf}
    hi = {"x": math.inf, "y": math.inf, "z": math.inf}
    inside, outside = [], []
    for token in re.findall(r"-?\d+", str(region)):
        surface, negative = surfaces[abs(int(token))], token.startswith("-")
        if isinstance(surface, openmc.ZCylinder):
            (inside if negative else outside).append(surface.r)
        else:
            axis = {openmc.XPlane: "x", openmc.YPlane: "y", openmc.ZPlane: "z"}[type(surface)]
            value = {"x": "x0", "y": "y0", "z": "z0"}[axis]
            if negative:
                hi[axis] = min(hi[axis], getattr(surface, value))
            else:
                lo[axis] = max(lo[axis], getattr(surface, value))
    if not (inside or outside) or all(math.isinf(lo[a]) and math.isinf(hi[a]) for a in "xy"):
        return None
    big = 1e6
    x0, x1 = max(lo["x"], -big), min(hi["x"], big)
    y0, y1 = max(lo["y"], -big), min(hi["y"], big)
    area = disk_rect_area(inside[0], x0, x1, y0, y1) if inside else (x1 - x0) * (y1 - y0)
    if outside:
        area -= disk_rect_area(outside[0], x0, x1, y0, y1)
    return area * (hi["z"] - lo["z"])


# ------------------------------------------------------------------------------
# Geometry plots and mesh tallies, following the MCNP plot command files ../mcnp/plot_<plane>_<view>
# ------------------------------------------------------------------------------
PLOTS = [(plane, view) for view in ["pit", "cavity"] for plane in ["xy", "xz"]]


@functools.cache
def read_plot(plane, view):
    """Return (origin, extent, {material number: color}) of an MCNP plot command file (its last
    origin and extent)."""
    origin = extent = None
    shades = {}
    for line in open(os.path.join(MCNP_DIR, f"plot_{plane}_{view}")):
        f = line.split("$")[0].split()
        if f and f[0] == "origin":
            origin = [float(v) for v in f[1:4]]
        elif f and f[0] == "extent":
            extent = float(f[1])
        elif f and f[0] == "shade":
            shades[int(f[1])] = f[2].lower()
    return origin, extent, shades


def rgb(name):
    """An MCNP plot color (an X11 color name, such as DodgerBlue) as 0-255 RGB."""
    return tuple(round(255 * v) for v in matplotlib.colors.to_rgb(name.lower()))


def get_plots(mats):
    plots = []
    for plane, view in PLOTS:
        origin, extent, shades = read_plot(plane, view)
        plot = openmc.SlicePlot()
        plot.basis = plane
        plot.origin = origin
        plot.width = (2 * extent, 2 * extent)
        plot.pixels = (1000, 1000)
        plot.color_by = "material"
        plot.colors = {mats[m]: rgb(color) for m, color in shades.items() if m in mats}
        plot.filename = os.path.join(IMAGES_DIR, f"artemis_openmc_geom_{plane}_{view}")
        plots.append(plot)
    return openmc.Plots(plots)


# Mesh tallies: over the window of each plot, 4 inches deep across the plot plane; 1-inch bins for
# the pit plots and half-inch bins for the cavity plots
MESH_DEPTH = 4 * 2.54
MESH_BINS = {"pit": 2.54, "cavity": 1.27}
# Energy groups of the neutron mesh tallies [eV]: below 0.5 eV, 0.5 eV to 0.1 MeV, 0.1 to 10 MeV,
# and above 10 MeV (up to the top of the cassette spectrum's bins), as emesh in the MCNP deck
NEUTRON_GROUPS = [0.0, 0.5, 1e5, 1e7, 1e9]


def get_tallies(geom):

    cells = geom.get_all_cells()
    cassettes = openmc.CellFilter([cells[106], cells[206]])
    # One filter per particle, shared by the tallies (OpenMC merges equal filters when it writes
    # them, and warns about it)
    only = {p: openmc.ParticleFilter(p) for p in ["neutron", "photon"]}
    groups = openmc.EnergyFilter(NEUTRON_GROUPS)
    tallies = []

    t = openmc.Tally(14, name="Neutron energy spectrum in the cassette regions (CCFE-709)")
    t.filters = [cassettes, only["neutron"], openmc.EnergyFilter.from_group_structure("CCFE-709")]
    t.scores = ["flux"]
    tallies.append(t)

    t = openmc.Tally(24, name="Photon energy spectrum in the cassette regions (VITAMIN-J-42)")
    t.filters = [
        cassettes,
        only["photon"],
        openmc.EnergyFilter.from_group_structure("VITAMIN-J-42"),
    ]
    t.scores = ["flux"]
    tallies.append(t)

    number = 34
    for particle in ["neutron", "photon"]:
        for plane, view in PLOTS:
            origin, extent, _ = read_plot(plane, view)
            mesh = openmc.RegularMesh(mesh_id=number)
            lower, upper, dims = [], [], []
            for axis, c0 in zip("xyz", origin):
                half, n = (
                    (extent, round(2 * extent / MESH_BINS[view]))
                    if axis in plane
                    else (MESH_DEPTH / 2, 1)
                )
                lower.append(c0 - half)
                upper.append(c0 + half)
                dims.append(n)
            mesh.lower_left, mesh.upper_right, mesh.dimension = lower, upper, dims
            t = openmc.Tally(
                number,
                name=f"{particle.capitalize()} flux in the {plane} plane (as plot_{plane}_{view})",
            )
            t.filters = [openmc.MeshFilter(mesh), only[particle]]
            if particle == "neutron":
                t.filters.append(groups)
            t.scores = ["flux"]
            tallies.append(t)
            number += 10

    return openmc.Tallies(tallies)


def get_settings(sources):
    # fmt: off
    settings                  = openmc.Settings()
    settings.run_mode         = 'fixed source'
    settings.source           = sources
    settings.photon_transport = True
    settings.batches          = 10
    settings.particles        = int(1e7)

    # A statepoint after every batch (each batch is 1e7 histories), with the source bank, so a run
    # that stops can be continued with --restart
    settings.statepoint       = {'batches': list(range(1, settings.batches + 1))}
    settings.sourcepoint      = {'write': True}

    # The tallies are in the statepoints: no tallies.out
    settings.output           = {'tallies': False}
    # fmt: on

    return settings


def statepoints():
    """Statepoint files in this folder, by batch."""
    files = glob.glob(os.path.join(HERE, "statepoint.*.h5"))
    return sorted(files, key=lambda f: int(f.split(".")[-2]))


def main():

    parser = argparse.ArgumentParser(description="Simplified OpenMC model of Artemis")
    parser.add_argument("--plot", action="store_true", help="make the geometry plots only")
    parser.add_argument("--inputs", action="store_true", help="write model.xml only")
    parser.add_argument(
        "--restart", action="store_true", help="continue from the latest statepoint"
    )
    parser.add_argument(
        "--threads", type=int, default=None, help="number of threads (default: all)"
    )
    args = parser.parse_args()

    check_source(SOURCE_H5)

    # fmt: off
    mats     = get_materials()
    geometry = get_geometry(mats)
    model    = openmc.Model(geometry=geometry, materials=openmc.Materials(mats.values()),
                            settings=get_settings(get_openmc_sources(SOURCE_H5)),
                            tallies=get_tallies(geometry), plots=get_plots(mats))
    # fmt: on

    if args.plot:
        model.plot_geometry(cwd=HERE)
        return
    if args.inputs:
        model.export_to_model_xml(os.path.join(HERE, "model.xml"))
        return

    restart = None
    if args.restart:
        if not statepoints():
            raise SystemExit("No statepoint to restart from")
        restart = statepoints()[-1]
    model.run(cwd=HERE, threads=args.threads, restart_file=restart)  # writes model.xml first

    # Keep only the final statepoint
    for f in statepoints()[:-1]:
        os.remove(f)


if __name__ == "__main__":
    main()
