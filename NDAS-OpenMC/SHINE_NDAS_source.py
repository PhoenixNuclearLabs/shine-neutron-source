# ******************************************************************************
# Copyright 2026, SHINE Technologies. All rights reserved.
# ******************************************************************************

import os
import h5py
import numpy as np
from itertools import product

import openmc

def get_sources(fname_h5, length=None, diameter=None, strength=None):

    # Load data from HDF5
    with h5py.File(fname_h5, 'r') as hf:
        x_grid      = hf['x_grid'     ][()]  # cm
        r_grid      = hf['radius_grid'][()]  # cm
        a_grid      = hf['angle_grid' ][()]  # degrees
        E_grid_DD   = hf['EN_grid_DD' ][()]  # MeV
        E_grid_DT   = hf['EN_grid_DT' ][()]  # MeV
        neutrons_DD = hf['neutrons_DD'][()]  # n/s
        neutrons_DT = hf['neutrons_DT'][()]  # n/s

    # Scale length
    if length is not None: x_grid = np.linspace(length / 2, -length / 2, x_grid.shape[0])

    # Account for smaller diameter if requested (cannot request larger diameter)
    # Note: beam scraping is not re-calculated here!
    if diameter is not None: r_grid[r_grid > diameter / 2] = diameter / 2

    # Scale neutron source to correct strength
    if strength is not None:
        factor = strength / (np.sum(neutrons_DD) + np.sum(neutrons_DT))
        neutrons_DD *= factor
        neutrons_DT *= factor

    # Calculate radius midpoints
    r_mids = np.sqrt((r_grid[:-1]**2 + r_grid[1:]**2) / 2)

    # Calculate direction cosines
    c_grid = np.cos(a_grid * np.pi / 180)

    # Neutrons_DD and Neutrons_DT have dimensions (position, angle, energy)
    num_x_ints = neutrons_DT.shape[0]
    num_a_ints = neutrons_DT.shape[1]

    E_grid   = [E_grid_DD  , E_grid_DT  ]
    neutrons = [neutrons_DD, neutrons_DT]

    sources = []

    # Create an IndependentSource for each position, angle, and fusion reaction
    for i, j, k in product(range(num_x_ints), range(num_a_ints), range(2)):

        total = np.sum(neutrons[k][i, j, :])
        if total == 0: continue

        source = openmc.IndependentSource()
        source.particle = 'neutron'
        source.strength = total

        # Spatial distribution - uniform within small cylindrical element
        source.space = openmc.stats.CylindricalIndependent(
            origin = (0, 0, 0),
            z   = openmc.stats.Tabular(x=x_grid[i:i + 2][::-1], p=[1], interpolation='histogram'),
            r   = openmc.stats.PowerLaw(a=0, b=r_mids[i], n=1),
            phi = openmc.stats.Uniform(a=0, b=2 * np.pi),
        )

        # Angular distribution - uniform within small solid angle
        source.angle = openmc.stats.PolarAzimuthal(
            reference_uvw = (0, 0, -1),
            mu  = openmc.stats.Tabular(x=c_grid[j:j + 2][::-1], p=[1], interpolation='histogram'),
            phi = openmc.stats.Uniform(a=0, b=2 * np.pi),
        )

        # Energy distribution - dependent on space and angle
        inds_nz = np.nonzero(neutrons[k][i, j, :])[0]
        source.energy = openmc.stats.Tabular(
            x = E_grid[k][inds_nz[0]:inds_nz[-1] + 2] * 1e6,
            p = np.append(neutrons[k][i, j, inds_nz], 0),
            interpolation = 'linear-linear',
        )

        sources.append(source)

    return sources
