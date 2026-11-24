#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies, LLC. MIT License: see ../../LICENSE.
# ******************************************************************************

# A very simplified model of Apollo, SHINE's neutron generator for UKAEA (the UK Atomic Energy
# Authority), in its bunker, the same as ../../mcnp/apollo/apollo.i: the same dimensions, cells,
# surfaces, materials, and tallies, with the same cell, surface, material, and tally numbers. It
# holds the target chamber with its double wall, hanging from the ceiling slab of the main
# irradiation room; the upstream section of the beamline (the beam tube runs up to the top of the
# slab) and the polyethylene around it in the slab's beamline penetration; the main and side
# irradiation rooms and the upper room (air); the concrete floor, walls, ceiling slab, and roof;
# and 3 ft of earth below the floor. The backstop is at the origin and the beam travels down (-z);
# the floor is at z = -308 cm and the ceiling slab runs from z = 146 to 266 cm.
#
# Differences from the MCNP model:
# - There is no cell for the void outside the model (MCNP cell 999): the outer surfaces of the
#   model are vacuum boundaries instead.
# - Tallies are already absolute rates: OpenMC scales them by the total strength of the sources,
#   which beam-target gives in n/s (their sum is the total neutron yield), as the MCNP model does
#   with factor. But flux tallies are not divided by volume: divide each mesh tally by the volume
#   of a mesh element. The cylindrical meshes of the peak tallies start at r = 4.9 cm, so each has
#   the one bin that is the second radial bin in MCNP.
#
# The neutron source comes from ../../source/Apollo/source.h5, made by beam-target
# (https://github.com/lucasjacobson/beam-target), and is read with
# ../../source/Apollo/openmc_source.py, which beam-target writes with it, so beam-target need not
# be installed to run this script. The geometry plots and mesh tallies follow the MCNP plot command
# files ../../mcnp/apollo/plot_*. Set OPENMC_CROSS_SECTIONS to the cross_sections.xml of an
# ENDF/B-VIII.1 library, the data the MCNP model uses.
#
# Run it from anywhere; OpenMC's files (model.xml, statepoints, ...) go in this folder:
#   python apollo.py               write model.xml and run (a statepoint after every batch)
#   python apollo.py --restart     continue the run from its latest statepoint
#   python apollo.py --plot        make the geometry plots (../../images/apollo/*_openmc_geom_*)
#   python apollo.py --inputs      only write model.xml
# Add --threads N to use N threads (default: all).

import argparse
import functools
import glob
import importlib.util
import math
import os
import re

import h5py
import matplotlib.colors
import numpy as np

import openmc

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE_DIR = os.path.join(HERE, "..", "..", "source", "Apollo")
SOURCE_H5 = os.path.join(SOURCE_DIR, "source.h5")


def get_openmc_sources(path):
    """The OpenMC sources of a source.h5 file, from the openmc_source.py that beam-target wrote with
    it (so beam-target need not be installed)."""
    spec = importlib.util.spec_from_file_location(
        "openmc_source", os.path.join(SOURCE_DIR, "openmc_source.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.get_openmc_sources(path)


MCNP_DIR = os.path.join(HERE, "..", "..", "mcnp", "apollo")
IMAGES_DIR = os.path.join(HERE, "..", "..", "images", "apollo")

# ------------------------------------------------------------------------------
# Dimensions [cm]
# ------------------------------------------------------------------------------
# fmt: off
# Target chamber and upstream section: radii
R_BEAM_IN, R_BEAM_OUT = 3.9624, 4.5974   # beam tube
R_BORE = 4.1783                          # target chamber inner radius (3.29 inches)
R_WALL_1, R_WATER, R_WALL_2 = 4.45807, 4.5974, 4.826
R_UP_IN, R_UP_OUT = 12.065, 12.7         # upstream tube
R_HOLE = 27.5                            # beamline penetration in the ceiling slab

# Target chamber and upstream section: heights
Z_END = -3.52552                         # bottom of the target chamber end plate
Z_END_TOP = -1.18872                     # top of the end plate
Z_PLATE = -0.9525                        # bottom of the steel plate below the beam stop
Z_STOP = -0.635                          # bottom of the beam stop
Z_WALLS = 129.54                         # top of the target chamber walls
Z_UP0, Z_UP1, Z_UP2 = 130.175, 208.026, 208.661  # upstream tube: bottom plate, top, top plate

# Bunker: outside of the walls, earth, floor, ceiling slab and its penetration, and roof
XB0, XB1, YB = -1067.6, 1412.4, 750.0
Z_EARTH, Z_FLOOR0, Z_FLOOR = round(-428.0 - 3 * 30.48, 5), -428.0, -308.0  # earth: 3 ft
Z_SLAB0, Z_HOLE, Z_FILL, Z_SLAB1 = 146.0, 162.0, 222.0, 266.0
Z_ROOF0, Z_ROOF1 = 1556.0, 1676.0
# Main irradiation room, in three boxes (x0, x1, y0, y1), the side room, and the upper room
ROOM_A = (-767.6, -93.4, -450.0, 450.0)
ROOM_B = (-93.4, 392.4, -450.0, 350.0)   # holds the target chamber
ROOM_C = (392.4, 592.4, -60.0, 350.0)
SIDE_ROOM = (792.4, 1292.4, -630.0, 630.0)
UPPER_ROOM = (-947.6, 1292.4, -630.0, 630.0)
# Polyethylene fill of the penetration (above the hole), in two boxes: x0, x1, and y1, and y0 of
# the lower and upper box
FILL_X0, FILL_X1, FILL_Y1 = -50.0, 105.0, 40.0
FILL_Y0_LOW, FILL_Y0_UP = -67.5, -89.5
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
# Polyethylene, Non-borated (C2H4), rho = 0.93 g/cm^3: here HDPE pellets at 75% of 0.95 g/cm^3
# (0.7125 g/cm^3) in the beamline penetration, as in the full model
COMPENDIUM[275] = ('Polyethylene, Non-borated (C2H4)', 0.7125, 'c_H_in_CH2', [
    ('H1', 6.6659000000e-01), ('H2', 7.6666666667e-05), ('C12', 3.2976666667e-01),
    ('C13', 3.5666666667e-03),
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
    rooms = [ROOM_A, ROOM_B, ROOM_C, SIDE_ROOM, UPPER_ROOM]
    room_x, room_y = [v for r in rooms for v in r[:2]], [v for r in rooms for v in r[2:]]
    groups = [
        ("cz", 101, [R_BEAM_IN, R_BORE, R_WALL_1, R_WATER, R_WALL_2, R_UP_IN, R_UP_OUT, R_HOLE]),
        ("px", 201, [XB0, XB1, FILL_X0, FILL_X1] + room_x),
        ("py", 301, [-YB, YB, FILL_Y0_LOW, FILL_Y0_UP, FILL_Y1] + room_y),
        (
            "pz",
            401,
            [
                Z_EARTH,
                Z_FLOOR0,
                Z_FLOOR,
                Z_END,
                Z_END_TOP,
                Z_PLATE,
                Z_STOP,
                0.0,
                Z_WALLS,
                Z_UP0,
                Z_SLAB0,
                Z_HOLE,
                Z_UP1,
                Z_UP2,
                Z_FILL,
                Z_SLAB1,
                Z_ROOF0,
                Z_ROOF1,
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
    outer = [("px", XB0), ("px", XB1), ("py", -YB), ("py", YB), ("pz", Z_EARTH), ("pz", Z_ROOF1)]
    for key in outer:
        surfaces[key].boundary_type = "vacuum"
    return surfaces


# Material numbers (None: void)
VOID, CU, SS, H2O, AIR, CONC, EARTH, POLY = None, 112, 331, 392, 4, 108, 116, 275

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
    (13, VOID, None,       R_BEAM_IN,  Z_UP1,     Z_SLAB1,   'Vacuum in beam tube'),
    (14, SS,   R_BEAM_IN,  R_BEAM_OUT, Z_UP1,     Z_SLAB1,   'Beam tube'),
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

    # Irradiation rooms
    xa0, xa1, ya0, ya1 = ROOM_A
    xb0, xb1, yb0, yb1 = ROOM_B
    xc0, xc1, yc0, yc1 = ROOM_C
    xs0, xs1, ys0, ys1 = SIDE_ROOM
    # fmt: off
    cell(101, AIR, box(xa0, xa1, ya0, ya1, Z_FLOOR, Z_SLAB0), 'Air in main irradiation room (-x)')
    cell(102, AIR, box(xb0, xb1, yb0, yb1, Z_FLOOR, Z_END),
         'Air in main irradiation room, below chamber')
    cell(103, AIR, box(xb0, xb1, yb0, yb1, Z_END, Z_WALLS) & +S['cz', R_WALL_2],
         'Air in main irradiation room, around chamber')
    cell(104, AIR, box(xb0, xb1, yb0, yb1, Z_WALLS, Z_SLAB0) & +S['cz', R_UP_OUT],
         'Air in main irradiation room, around upstream tube')
    cell(105, AIR, box(xc0, xc1, yc0, yc1, Z_FLOOR, Z_SLAB0), 'Air in main irradiation room (+x)')
    cell(106, AIR, box(xs0, xs1, ys0, ys1, Z_FLOOR, Z_SLAB0), 'Air in side irradiation room')

    # Beamline penetration in the ceiling slab: the hole, and the polyethylene above it
    cell(201, AIR, tube(R_UP_OUT, R_HOLE, Z_SLAB0, Z_HOLE), 'Air in penetration hole')
    fill = (FILL_X0, FILL_X1)
    cell(202, POLY, box(*fill, FILL_Y0_LOW, FILL_Y1, Z_HOLE, Z_UP2) & +S['cz', R_UP_OUT],
         'Penetration polyethylene, around upstream tube')
    cell(203, POLY, box(*fill, FILL_Y0_LOW, FILL_Y1, Z_UP2, Z_FILL) & +S['cz', R_BEAM_OUT],
         'Penetration polyethylene, around beam tube (lower)')
    cell(204, POLY, box(*fill, FILL_Y0_UP, FILL_Y1, Z_FILL, Z_SLAB1) & +S['cz', R_BEAM_OUT],
         'Penetration polyethylene, around beam tube (upper)')

    # Upper room
    xu0, xu1, yu0, yu1 = UPPER_ROOM
    cell(301, AIR, box(xu0, xu1, yu0, yu1, Z_SLAB1, Z_ROOF0), 'Air in upper room')

    # Concrete: floor, walls, ceiling slab, upper walls, and roof
    z = (Z_FLOOR, Z_SLAB0)
    cell(401, CONC, box(XB0, XB1, -YB, YB, Z_FLOOR0, Z_FLOOR), 'Floor')
    cell(402, CONC, box(XB0, xa0, -YB, YB, *z), 'Walls, -x end')
    cell(403, CONC, box(xa0, xa1, -YB, ya0, *z), 'Walls, beside main room (-x), -y side')
    cell(404, CONC, box(xa0, xa1, ya1, YB, *z), 'Walls, beside main room (-x), +y side')
    cell(405, CONC, box(xb0, xb1, -YB, yb0, *z), 'Walls, beside main room, -y side')
    cell(406, CONC, box(xb0, xb1, yb1, YB, *z), 'Walls, beside main room, +y side')
    cell(407, CONC, box(xc0, xc1, -YB, yc0, *z), 'Walls, beside main room (+x), -y side')
    cell(408, CONC, box(xc0, xc1, yc1, YB, *z), 'Walls, beside main room (+x), +y side')
    cell(409, CONC, box(xc1, xs0, -YB, YB, *z), 'Walls, between main and side rooms')
    cell(410, CONC, box(xs0, xs1, -YB, ys0, *z), 'Walls, beside side room, -y side')
    cell(411, CONC, box(xs0, xs1, ys1, YB, *z), 'Walls, beside side room, +y side')
    cell(412, CONC, box(xs1, XB1, -YB, YB, *z), 'Walls, +x end')
    cell(421, CONC, box(XB0, XB1, -YB, YB, Z_SLAB0, Z_HOLE) & +S['cz', R_HOLE],
         'Ceiling slab, around penetration hole')
    for n, z0, z1, y0, part in [(422, Z_HOLE, Z_FILL, FILL_Y0_LOW, 'lower'),
                                (426, Z_FILL, Z_SLAB1, FILL_Y0_UP, 'upper')]:
        name = f'Ceiling slab, {{}} of penetration polyethylene ({part})'
        cell(n,     CONC, box(XB0, FILL_X0, -YB, YB, z0, z1), name.format('-x'))
        cell(n + 1, CONC, box(FILL_X1, XB1, -YB, YB, z0, z1), name.format('+x'))
        cell(n + 2, CONC, box(*fill, -YB, y0, z0, z1), name.format('-y'))
        cell(n + 3, CONC, box(*fill, FILL_Y1, YB, z0, z1), name.format('+y'))
    cell(431, CONC, box(XB0, xu0, -YB, YB, Z_SLAB1, Z_ROOF0), 'Upper walls, -x end')
    cell(432, CONC, box(xu0, xu1, -YB, yu0, Z_SLAB1, Z_ROOF0), 'Upper walls, -y side')
    cell(433, CONC, box(xu0, xu1, yu1, YB, Z_SLAB1, Z_ROOF0), 'Upper walls, +y side')
    cell(434, CONC, box(xu1, XB1, -YB, YB, Z_SLAB1, Z_ROOF0), 'Upper walls, +x end')
    cell(441, CONC, box(XB0, XB1, -YB, YB, Z_ROOF0, Z_ROOF1), 'Roof')

    # Earth
    cell(501, EARTH, box(XB0, XB1, -YB, YB, Z_EARTH, Z_FLOOR0), 'Earth below floor')
    # fmt: on

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
def disk_rect_area(r, x0, x1, y0, y1, xc=0.0):
    """Area of the rectangle [x0, x1] x [y0, y1] inside the disk of radius r about (xc, 0)."""
    x0, x1 = max(x0 - xc, -r), min(x1 - xc, r)
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
            (inside if negative else outside).append((surface.r, surface.x0))
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
    if inside:
        area = disk_rect_area(inside[0][0], x0, x1, y0, y1, inside[0][1])
    else:
        area = (x1 - x0) * (y1 - y0)
    if outside:
        area -= disk_rect_area(outside[0][0], x0, x1, y0, y1, outside[0][1])
    return area * (hi["z"] - lo["z"])


# ------------------------------------------------------------------------------
# Geometry plots and mesh tallies, following the MCNP plot command files
# ../../mcnp/apollo/plot_<plane>_<view>
# ------------------------------------------------------------------------------
PLOTS = [(plane, view) for view in ["room", "chamber"] for plane in ["xy", "xz"]]


@functools.cache
def read_plot(plane, view):
    """Return (origin, (horizontal, vertical) half-widths, {material number: color}) of an MCNP plot
    command file (its last origin and extent; one extent value is a square window)."""
    origin = extent = None
    shades = {}
    for line in open(os.path.join(MCNP_DIR, f"plot_{plane}_{view}")):
        f = line.split("$")[0].split()
        if f and f[0] == "origin":
            origin = [float(v) for v in f[1:4]]
        elif f and f[0] == "extent":
            extent = [float(v) for v in f[1:3]]
            extent = (extent[0], extent[-1])
        elif f and f[0] == "shade":
            shades[int(f[1])] = f[2].lower()
    return origin, extent, shades


def rgb(name):
    """An MCNP plot color (an X11 color name, such as DodgerBlue) as 0-255 RGB."""
    return tuple(round(255 * v) for v in matplotlib.colors.to_rgb(name.lower()))


def get_plots(mats):
    plots = []
    for plane, view in PLOTS:
        origin, (eh, ev), shades = read_plot(plane, view)
        plot = openmc.SlicePlot()
        plot.basis = plane
        plot.origin = origin
        plot.width = (2 * eh, 2 * ev)
        plot.pixels = (1000, round(1000 * ev / eh))  # square pixels
        plot.color_by = "material"
        plot.colors = {mats[m]: rgb(color) for m, color in shades.items() if m in mats}
        plot.filename = os.path.join(IMAGES_DIR, f"apollo_openmc_geom_{plane}_{view}")
        plots.append(plot)
    return openmc.Plots(plots)


# Mesh tallies: over the window of each plot, 4 inches deep across the plot plane; 2-inch bins for
# the room plots and half-inch bins for the chamber plots
MESH_DEPTH = 4 * 2.54
MESH_BINS = {"room": 5.08, "chamber": 1.27}
# Energy groups of the neutron mesh tallies [eV]: below 0.5 eV, 0.5 eV to 0.1 MeV, 0.1 to 10 MeV,
# and above 10 MeV (up to the top of the spectrum's bins), as emesh in the MCNP deck
NEUTRON_GROUPS = [0.0, 0.5, 1e5, 1e7, 1e9]
# Peak values: a cylindrical shell against the chamber's outer wall, from r = 4.9 to 5 cm, 5 cm
# tall, centered at the height of the peak flux (the height of the plot_xy planes)
R_PEAK, Z_PEAK, H_PEAK = (4.9, 5.0), 55.88, 5.0
# Responses in silicon, from the tables in data/ (from MCNP's own ACE files for silicon, so they are
# the same data as the MCNP deck's fm cards with m311 use):
# - 1-MeV(Si) equivalent neutron flux: the flux times the damage energy cross section (MT 444) over
#   the 1 MeV reference of 95 MeV-mb
# - neutron and photon kerma in silicon [rad(Si)/s]: the flux times the heating cross section (the
#   total cross section times the heating number) times RAD
DATA_DIR = os.path.join(HERE, "..", "..", "data")
D_REF = 95e3  # [eV-b]
A_SI = 28.0855  # [g/mol]
RAD = 6.02214076e23 / A_SI * 1e-24 * 1.602176634e-8 * 1e-6  # [rad-cm^2 per eV-b]


@functools.cache
def silicon_filters():
    """EnergyFunctionFilters of the responses in silicon: {"si", "nkerma", "pkerma"}."""
    neutron = np.loadtxt(os.path.join(DATA_DIR, "silicon_neutron.txt"))
    photon = np.loadtxt(os.path.join(DATA_DIR, "silicon_photon.txt"))
    # The photon data give each photoelectric edge twice (the values below and above the jump);
    # OpenMC needs increasing energies, so the second of each pair moves up by a relative 1e-9
    for table in [neutron, photon]:
        repeat = np.concatenate([[False], np.diff(table[:, 0]) == 0])
        table[repeat, 0] *= 1 + 1e-9
        assert np.all(np.diff(table[:, 0]) > 0)
    return {
        "si": openmc.EnergyFunctionFilter(neutron[:, 0], neutron[:, 1] / D_REF),
        "nkerma": openmc.EnergyFunctionFilter(neutron[:, 0], neutron[:, 2] * RAD),
        "pkerma": openmc.EnergyFunctionFilter(
            photon[:, 0], photon[:, 1] * RAD, interpolation="log-log"
        ),
    }


# The quantities tallied: name, particle, and the response filter (None for the flux)
QUANTITIES = {
    "neutron": ("Neutron flux", "neutron"),
    "photon": ("Photon flux", "photon"),
    "si": ("1-MeV(Si) equivalent neutron flux", "neutron"),
    "nkerma": ("Neutron kerma in silicon", "neutron"),
    "pkerma": ("Photon kerma in silicon", "photon"),
}
# Tally numbers, as in the MCNP deck (see its header). The tens digit q of the peak tallies (2q4) is
# the quantity. The maps take a hundred per quantity (4xx to 8xx), with the window in the tens
# digit: 1 and 2 for the xy and xz planes of the room view, 4 and 5 for those of the chamber view
# (3 and 6 are kept for yz planes).
Q = {"neutron": 1, "nspectrum": 2, "pspectrum": 3, "si": 4, "nkerma": 5, "pkerma": 6}
MAP_BASE = {"neutron": 400, "photon": 500, "si": 600, "nkerma": 700, "pkerma": 800}
WINDOW = {"xy": 1, "xz": 2, "yz": 3}
ZOOM = {"room": 0, "chamber": 3}
# Mesh tallies: each quantity in each plot window
MAPS = ["neutron", "photon", "si", "nkerma", "pkerma"]


def shell(mesh_id, r, z0, z1):
    """A one-bin cylindrical mesh filter over the shell r[0] < r < r[1], z0 < z < z1."""
    mesh = openmc.CylindricalMesh(
        r_grid=list(r), z_grid=[z0, z1], phi_grid=[0.0, 2 * math.pi], mesh_id=mesh_id
    )
    return openmc.MeshFilter(mesh)


def get_tallies():

    # One filter per particle, shared by the tallies (OpenMC merges equal filters when it writes
    # them, and warns about it)
    only = {p: openmc.ParticleFilter(p) for p in ["neutron", "photon"]}
    groups = openmc.EnergyFilter(NEUTRON_GROUPS)
    ccfe = openmc.EnergyFilter.from_group_structure("CCFE-709")
    vitamin_j = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")
    # The filter of each quantity's energy dependence
    energy = {"neutron": groups, "photon": None, **silicon_filters()}
    tallies = []

    def tally(number, name, *filters):
        t = openmc.Tally(number, name=name)
        t.filters = [f for f in filters if f is not None]
        t.scores = ["flux"]
        tallies.append(t)

    # The peak (2q4): one-bin cylindrical meshes, for neutrons in the four groups, the neutron
    # spectrum, the photon spectrum, the 1-MeV(Si) equivalent neutron flux, and the neutron and
    # photon kerma in silicon
    # The mesh is numbered as its first tally
    peak = shell(214, R_PEAK, Z_PEAK - H_PEAK / 2, Z_PEAK + H_PEAK / 2)
    for q, what, particle, energy_filter in [
        ("neutron", "neutron flux", "neutron", groups),
        ("nspectrum", "neutron energy spectrum (CCFE-709)", "neutron", ccfe),
        ("pspectrum", "photon energy spectrum (VITAMIN-J-42)", "photon", vitamin_j),
        ("si", "1-MeV(Si) equivalent neutron flux", "neutron", energy["si"]),
        ("nkerma", "neutron kerma in silicon", "neutron", energy["nkerma"]),
        ("pkerma", "photon kerma in silicon", "photon", energy["pkerma"]),
    ]:
        tally(200 + 10 * Q[q] + 4, f"Peak {what} around the chamber", peak, only[particle],
              energy_filter)  # fmt: skip

    for quantity in MAPS:
        name, particle = QUANTITIES[quantity]
        for plane, view in PLOTS:
            number = MAP_BASE[quantity] + 10 * (ZOOM[view] + WINDOW[plane]) + 4
            origin, extent, _ = read_plot(plane, view)
            mesh = openmc.RegularMesh(mesh_id=number)
            lower, upper, dims = [], [], []
            for axis, c0 in zip("xyz", origin):
                if axis in plane:
                    half = extent[plane.index(axis)]
                    n = round(2 * half / MESH_BINS[view])
                else:
                    half, n = MESH_DEPTH / 2, 1
                lower.append(c0 - half)
                upper.append(c0 + half)
                dims.append(n)
            mesh.lower_left, mesh.upper_right, mesh.dimension = lower, upper, dims
            tally(number, f"{name} in the {plane} plane (as plot_{plane}_{view})",
                  openmc.MeshFilter(mesh), only[particle], energy[quantity])  # fmt: skip

    return openmc.Tallies(tallies)


def get_settings(sources):
    # fmt: off
    settings                  = openmc.Settings()
    settings.run_mode         = 'fixed source'
    settings.source           = sources
    settings.photon_transport = True
    settings.batches          = 10
    settings.particles        = int(1e6)

    # A statepoint after every batch (each batch is 1e6 histories), with the source bank, so a run
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

    parser = argparse.ArgumentParser(description="Simplified OpenMC model of Apollo")
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
                            tallies=get_tallies())
    # fmt: on

    # The geometry plots only with --plot: OpenMC checks that their folder (../../images/apollo/)
    # exists before any run that holds them, and doesn't make it
    if args.plot:
        os.makedirs(IMAGES_DIR, exist_ok=True)
        model.plots = get_plots(mats)
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
