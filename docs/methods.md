# Methods

What the models of every system have in common: how their tallies are numbered, how the responses
in silicon are calculated, how the results are normalized, how MCNP and OpenMC compare, and the
limitations of the models.

[Back to the README](../README.md) · [Artemis](artemis.md) · [Apollo](apollo.md)

## Contents

- [Tallies](#tallies)
- [Responses in silicon](#responses-in-silicon)
- [Units and normalization](#units-and-normalization)
- [MCNP and OpenMC](#mcnp-and-openmc)
- [Limitations](#limitations)

## Tallies

Every tally is a track-length (F4) or mesh tally, and a tally number means the same thing in every
system's models, in both codes; a model leaves out those it has no use for. Each number ends in 4:

| Tallies | What |
|---|---|
| 1q4 | Cell tallies in Artemis's cassette regions (Artemis only) |
| 2q4, 3q4 | Peak values: on Artemis's plates in its −x and +x cassette regions, or in Apollo's shell against the chamber wall (2q4 only) |
| 4w4 to 8w4 | Maps over the windows of the geometry plots: the neutron flux in four energy groups (4w4), the photon flux (5w4), the 1-MeV(Si) equivalent flux (6w4), and the neutron and photon kerma in silicon (7w4, 8w4) |

The tens digit q is the quantity: 1 the neutron flux in four groups, 2 the neutron spectrum
(CCFE-709 groups), 3 the photon spectrum (VITAMIN-J-42), 4 the 1-MeV(Si) equivalent flux, and 5
and 6 the neutron and photon kerma in silicon. The tens digit w of the maps is the plot window: 1
and 2 for the xy and xz planes of the whole view (Artemis's pit, Apollo's room), and 4 and 5 for
those of the zoomed view (Artemis's cavities, Apollo's chamber); 3 and 6 are kept for yz planes. The neutron energy groups are
below 0.5 eV, 0.5 eV to 0.1 MeV, 0.1 to 10 MeV, and above 10 MeV.

## Responses in silicon

- **1-MeV(Si) equivalent neutron flux:** the neutron flux weighted by silicon's displacement
  damage, the damage energy cross section (MT 444) of ENDF/B-VIII.1 silicon, over the 1 MeV
  reference of 95 MeV-mb. This is the ENDF/B-VIII.1 damage function, not the ASTM E722 table,
  which is based on older data.
- **Kerma in silicon:** the energy given to charged particles in a small piece of silicon, from
  the heating cross sections of ENDF/B-VIII.1 silicon (neutrons, from Lib81) and of eprdata25
  (photons), in rad(Si)/s. It equals the absorbed dose where there is charged-particle
  equilibrium, as in a device within a few millimeters of material.
- **Hardness factor:** the 1-MeV(Si) equivalent flux over the total neutron flux.

MCNP tallies these with `fm` cards on a silicon material that is in no cell. OpenMC tallies them
with energy-function filters from [`data/silicon_*.txt`](../data/), which hold the same data, taken
from the same MCNP ACE files, so the two codes use the same responses.

## Units and normalization

Fluxes are per cm² per second and kerma is in rad(Si) per second, at each system's full neutron
yield: 1.748845 × 10¹³ n/s for Artemis and 3.113519 × 10¹³ n/s for Apollo (`yield_total` in
`source/<System>/summary.toml`).

- **MCNP** tallies are per source neutron; the models multiply them by the yield with `fm` cards
  (cell tallies) and the `factor` keyword (mesh tallies). If a source changes, change the yield
  in these, including the constants of Artemis's silicon cell tallies 144, 154, and 164, which
  include it.
- **OpenMC** tallies are scaled by the source strength, which beam-target gives in n/s, so they
  are already absolute. Flux tallies are not divided by volume: `plot_results.py` divides the cell
  tallies by the cells' volumes and the mesh tallies by a mesh element's volume.

## MCNP and OpenMC

Each system's two models are the same, number for number: the same surfaces, cells, materials,
thermal scattering data, and tallies, with the same nuclear data (ENDF/B-VIII.1 neutrons and
eprdata25 photons) and the same responses in silicon.

- **Neutrons** agree within 0.6% in both systems at the peak, and the maps agree bin by bin within
  a few percent.
- **Photons** do not: OpenMC's photon flux and photon kerma are higher, by 2% to 4% near Apollo's
  chamber, 8% at Artemis's peak and 11% over its cassettes (photon kerma), and up to 16% across
  Apollo's room.
  The difference is mostly in the photon lines from oxygen-16 at 6 to 7.5 MeV, which neutrons above
  6 MeV excite in water, air, and concrete, and to a lesser extent below 0.1 MeV. Since both codes
  use the same evaluations, this points to how the photon production data were processed for each
  code; it is not yet resolved.

## Limitations

- **Simplified models.** Each leaves out the accelerator, the internal structure of the beam stop,
  and any parts or fixtures where parts go, which is left empty; parts there change the field
  around them.
  See each system's page for what else is left out.
- **The MCNP sources** cannot hold the dependence of energy and angle on position: an SDEF
  distribution can depend on only one other variable, so in `source.sdef` the neutron energy
  depends on the emission angle but not on the position along the chamber. The OpenMC sources
  (`source.h5`) keep the full correlation between position, angle, and energy.
- **The photon results** of the two codes differ (see [MCNP and OpenMC](#mcnp-and-openmc)).
- **Source assumptions** are beam-target's (straggling neglected, a flat beam profile, and
  others); see its README.
