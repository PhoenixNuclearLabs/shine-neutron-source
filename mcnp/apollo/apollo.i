SHINE Apollo neutron generator, simplified model
c
c  Author:       Lucas Jacobson
c  Created on:   November 23, 2026
c  Last updated: November 23, 2026
c
c  Copyright 2026, SHINE Technologies, LLC. MIT License: see ../../LICENSE.
c
c  A very simplified model of Apollo, SHINE's neutron generator for UKAEA
c  (the UK Atomic Energy Authority), in its bunker: the target chamber with
c  its double wall, hanging from the ceiling slab of the main irradiation room;
c  the upstream section of the beamline (the beam tube runs up to the top of
c  the slab) and the polyethylene around it in the slab's beamline
c  penetration; the main and side irradiation rooms and the upper room (air);
c  the concrete floor, walls, ceiling slab, and roof; and 3 ft of earth below
c  the floor. The irradiation rooms are empty. Every cell is a plain
c  intersection of half-spaces.
c
c  The backstop is at the origin and the beam travels down (-z) from the
c  top of the 137 cm source region. The floor is at z = -308 cm and the
c  ceiling slab runs from z = 146 to 266 cm.
c  Tallies are scaled by the total neutron yield (yield_total in
c  ../../source/Apollo/summary.toml). Responses in silicon use fm cards with
c  silicon (m311, which is in no cell): the 1-MeV(Si) equivalent neutron flux
c  (the damage energy cross section, MT 444, over the 1 MeV reference of
c  95 MeV-mb) and the neutron and photon kerma in silicon [rad(Si)/s]. The
c  OpenMC model uses the same data, from ../../data/silicon_*.txt.
c
c  Geometry plots: mcnp6 ip i=apollo.i com=plot_xy_room (or plot_xz_room,
c  plot_yz_room, or the zoomed plot_xy_chamber, plot_xz_chamber, and
c  plot_yz_chamber).
c
c  Tally numbers mean the same in every deck of this repository
c  (mcnp/<system>/<system>.i); a deck leaves out those it has no use for:
c    1q4         cell tallies in Artemis's cassette regions (Artemis only)
c    2q4         peak values: in Apollo, around the target chamber, in a
c                cylindrical shell from r = 4.9 to 5 cm (against the
c                chamber's outer wall), 5 cm tall, at the height of the peak
c                flux (the second radial bin; the first is inside the shell)
c    3q4         peak values in Artemis's +x cassette region (Artemis only)
c    4w4 to 8w4  maps over the windows of the geometry plots: the neutron
c                flux in four groups (4w4), the photon flux (5w4), the
c                1-MeV(Si) equivalent flux (6w4), and the neutron and photon
c                kerma in silicon (7w4, 8w4)
c  where q is the quantity (1 neutrons in four groups, 2 the neutron
c  spectrum, 3 the photon spectrum, 4 the 1-MeV(Si) equivalent flux, 5 and 6
c  the neutron and photon kerma in silicon) and w the window (1 to 3: the
c  xy, xz, and yz planes of the room view; 4 to 6: the chamber view).
c
c ******************************************************************************
c
c                                  CELL CARDS
c
c ******************************************************************************
c
c ------------------------------------------------------------------------------
c     Target chamber and upstream section
c ------------------------------------------------------------------------------
c
    1    0                  -102  408 -410                                  $ Vacuum in target chamber
    2  112  -8.96           -102  407 -408                                  $ Copper beam stop
    3  331  -8.03           -103  406 -407                                  $ Steel plate below beam stop
    4  331  -8.03       102 -103  407 -409                                  $ Target chamber inner wall
    5  392  -0.997      103 -104  406 -409                                  $ Cooling water jacket
    6  392  -0.997          -104  405 -406                                  $ Cooling water behind beam stop
    7  331  -8.03       104 -105  405 -409                                  $ Target chamber outer wall
    8  331  -8.03           -105  404 -405                                  $ Target chamber end plate
    9  331  -8.03       102 -107  409 -410                                  $ Upstream tube bottom plate
   10    0                  -106  410 -413                                  $ Vacuum in upstream tube
   11  331  -8.03       106 -107  410 -413                                  $ Upstream tube
   12  331  -8.03       104 -107  413 -414                                  $ Upstream tube top plate
   13    0                  -101  413 -416                                  $ Vacuum in beam tube
   14  331  -8.03       101 -104  413 -416                                  $ Beam tube
c
c ------------------------------------------------------------------------------
c     Irradiation rooms
c ------------------------------------------------------------------------------
c
  101    4  -0.001205   203 -204  303 -309  403 -411                        $ Air in main irradiation room (-x)
  102    4  -0.001205   204 -207  303 -308  403 -404                        $ Air in main irradiation room, below chamber
  103    4  -0.001205   204 -207  303 -308  404 -409  105  vol=51704847.48  $ Air in main irradiation room, around chamber
  104    4  -0.001205   204 -207  303 -308  409 -411  107  vol=6388673.99   $ Air in main irradiation room, around upstream tube
  105    4  -0.001205   207 -208  306 -308  403 -411                        $ Air in main irradiation room (+x)
  106    4  -0.001205   209 -210  302 -310  403 -411                        $ Air in side irradiation room
c
c ------------------------------------------------------------------------------
c     Beamline penetration in the ceiling slab
c ------------------------------------------------------------------------------
c
  201    4  -0.001205   107 -108  411 -412                                  $ Air in penetration hole
  202  275  -0.7125     205 -206  305 -307  412 -414  107  vol=753845.43    $ Penetration polyethylene, around upstream tube
  203  275  -0.7125     205 -206  305 -307  414 -415  104  vol=221375.36    $ Penetration polyethylene, around beam tube (lower)
  204  275  -0.7125     205 -206  304 -307  415 -416  104  vol=880268.36    $ Penetration polyethylene, around beam tube (upper)
c
c ------------------------------------------------------------------------------
c     Upper room
c ------------------------------------------------------------------------------
c
  301    4  -0.001205   202 -210  302 -310  416 -417                        $ Air in upper room
c
c ------------------------------------------------------------------------------
c     Concrete
c ------------------------------------------------------------------------------
c
  401  108  -2.3        201 -211  301 -311  402 -403                        $ Floor
  402  108  -2.3        201 -203  301 -311  403 -411                        $ Walls, -x end
  403  108  -2.3        203 -204  301 -303  403 -411                        $ Walls, beside main room (-x), -y side
  404  108  -2.3        203 -204  309 -311  403 -411                        $ Walls, beside main room (-x), +y side
  405  108  -2.3        204 -207  301 -303  403 -411                        $ Walls, beside main room, -y side
  406  108  -2.3        204 -207  308 -311  403 -411                        $ Walls, beside main room, +y side
  407  108  -2.3        207 -208  301 -306  403 -411                        $ Walls, beside main room (+x), -y side
  408  108  -2.3        207 -208  308 -311  403 -411                        $ Walls, beside main room (+x), +y side
  409  108  -2.3        208 -209  301 -311  403 -411                        $ Walls, between main and side rooms
  410  108  -2.3        209 -210  301 -302  403 -411                        $ Walls, beside side room, -y side
  411  108  -2.3        209 -210  310 -311  403 -411                        $ Walls, beside side room, +y side
  412  108  -2.3        210 -211  301 -311  403 -411                        $ Walls, +x end
  421  108  -2.3        201 -211  301 -311  411 -412  108  vol=59481986.73  $ Ceiling slab, around penetration hole
  422  108  -2.3        201 -205  301 -311  412 -415                        $ Ceiling slab, -x of penetration polyethylene (lower)
  423  108  -2.3        206 -211  301 -311  412 -415                        $ Ceiling slab, +x of penetration polyethylene (lower)
  424  108  -2.3        205 -206  301 -305  412 -415                        $ Ceiling slab, -y of penetration polyethylene (lower)
  425  108  -2.3        205 -206  307 -311  412 -415                        $ Ceiling slab, +y of penetration polyethylene (lower)
  426  108  -2.3        201 -205  301 -311  415 -416                        $ Ceiling slab, -x of penetration polyethylene (upper)
  427  108  -2.3        206 -211  301 -311  415 -416                        $ Ceiling slab, +x of penetration polyethylene (upper)
  428  108  -2.3        205 -206  301 -304  415 -416                        $ Ceiling slab, -y of penetration polyethylene (upper)
  429  108  -2.3        205 -206  307 -311  415 -416                        $ Ceiling slab, +y of penetration polyethylene (upper)
  431  108  -2.3        201 -202  301 -311  416 -417                        $ Upper walls, -x end
  432  108  -2.3        202 -210  301 -302  416 -417                        $ Upper walls, -y side
  433  108  -2.3        202 -210  310 -311  416 -417                        $ Upper walls, +y side
  434  108  -2.3        210 -211  301 -311  416 -417                        $ Upper walls, +x end
  441  108  -2.3        201 -211  301 -311  417 -418                        $ Roof
c
c ------------------------------------------------------------------------------
c     Earth
c ------------------------------------------------------------------------------
c
  501  116  -1.52       201 -211  301 -311  401 -402                        $ Earth below floor
c
c ------------------------------------------------------------------------------
c     Outside the model
c ------------------------------------------------------------------------------
c
  999    0             (-201:211:-301:311:-401:418)                         $ Rest of universe

c ******************************************************************************
c
c                                SURFACE CARDS
c
c ******************************************************************************
c
c ------------------------------------------------------------------------------
c     Cylinders along z
c ------------------------------------------------------------------------------
c
  101  cz       3.96240  $ Beam tube
  102  cz       4.17830  $ Target chamber inner radius
  103  cz       4.45807  $ Target chamber inner wall
  104  cz       4.59740  $ Cooling water jacket and beam tube
  105  cz       4.82600  $ Target chamber outer wall
  106  cz      12.06500  $ Upstream tube
  107  cz      12.70000  $ Upstream tube
  108  cz      27.50000  $ Beamline penetration
c
c ------------------------------------------------------------------------------
c     Planes normal to x
c ------------------------------------------------------------------------------
c
  201  px   -1067.60000  $ Bunker walls
  202  px    -947.60000  $ Upper room
  203  px    -767.60000  $ Main irradiation room
  204  px     -93.40000  $ Main irradiation room
  205  px     -50.00000  $ Penetration polyethylene
  206  px     105.00000  $ Penetration polyethylene
  207  px     392.40000  $ Main irradiation room
  208  px     592.40000  $ Main irradiation room
  209  px     792.40000  $ Side irradiation room
  210  px    1292.40000  $ Side irradiation room and upper room
  211  px    1412.40000  $ Bunker walls
c
c ------------------------------------------------------------------------------
c     Planes normal to y
c ------------------------------------------------------------------------------
c
  301  py    -750.00000  $ Bunker walls
  302  py    -630.00000  $ Side irradiation room and upper room
  303  py    -450.00000  $ Main irradiation room
  304  py     -89.50000  $ Penetration polyethylene (upper)
  305  py     -67.50000  $ Penetration polyethylene (lower)
  306  py     -60.00000  $ Main irradiation room
  307  py      40.00000  $ Penetration polyethylene
  308  py     350.00000  $ Main irradiation room
  309  py     450.00000  $ Main irradiation room
  310  py     630.00000  $ Side irradiation room and upper room
  311  py     750.00000  $ Bunker walls
c
c ------------------------------------------------------------------------------
c     Planes normal to z
c ------------------------------------------------------------------------------
c
  401  pz    -519.44000  $ Bottom of earth
  402  pz    -428.00000  $ Bottom of floor
  403  pz    -308.00000  $ Floor
  404  pz      -3.52552  $ Bottom of target chamber
  405  pz      -1.18872  $ Top of end plate (bottom of cooling water)
  406  pz      -0.95250  $ Bottom of steel plate below beam stop
  407  pz      -0.63500  $ Bottom of beam stop
  408  pz       0.00000  $ Backstop (top of beam stop)
  409  pz     129.54000  $ Top of target chamber walls
  410  pz     130.17500  $ Top of upstream tube bottom plate
  411  pz     146.00000  $ Ceiling (bottom of ceiling slab)
  412  pz     162.00000  $ Top of penetration hole
  413  pz     208.02600  $ Top of upstream tube
  414  pz     208.66100  $ Top of upstream tube top plate
  415  pz     222.00000  $ Penetration polyethylene (upper)
  416  pz     266.00000  $ Top of ceiling slab (end of beam tube)
  417  pz    1556.00000  $ Bottom of roof
  418  pz    1676.00000  $ Top of roof

c ******************************************************************************
c
c                                  DATA CARDS
c
c ******************************************************************************
c
c ------------------------------------------------------------------------------
c     Simulation parameters
c ------------------------------------------------------------------------------
c
  nps         1e7
  prdmp       j 1e7 1 2 j
c
  mode        n p
  imp:n       1 51r 0
  imp:p       1 51r 0
c
c ------------------------------------------------------------------------------
c     Source definition
c ------------------------------------------------------------------------------
c
  read        file = ../../source/Apollo/source.sdef
c
c ------------------------------------------------------------------------------
c     Tallies
c ------------------------------------------------------------------------------
c
c     Responses in silicon (m311, in no cell), with fm cards: the 1-MeV(Si) equivalent neutron
c     flux (the damage energy cross section, MT 444 [MeV-b], over the 1 MeV reference of 95 MeV-mb),
c     and the neutron and photon kerma [rad(Si)/s] (the total cross section times the heating
c     number [MeV-b], times 3.43541e-10 rad-cm^2 per MeV-b). The yield is in factor.
c
c ------------------------------------------------------------------------------
c     Mesh tallies
c ------------------------------------------------------------------------------
c
c     Peak values around the target chamber: a cylindrical shell from r = 4.9 to 5 cm (against
c     the chamber's outer wall), 5 cm tall, centered at z = 55.88 cm, the second radial bin
c     (the first is inside it): neutrons in the four groups, the neutron spectrum (CCFE-709),
c     the photon spectrum (VITAMIN-J-42), the 1-MeV(Si) equivalent neutron flux, and the
c     neutron and photon kerma in silicon
c
  fc214       Peak neutron flux around the chamber
  fmesh214:n  geom   = cyl
              origin = 0 0 53.38
              axs    = 0 0 1
              vec    = 1 0 0
              imesh  = 4.9 5     iints = 1 1
              jmesh  = 5         jints = 1
              kmesh  = 1         kints = 1
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
  fc224       Peak neutron spectrum around the chamber (CCFE-709)
  fmesh224:n  geom   = cyl
              origin = 0 0 53.38
              axs    = 0 0 1
              vec    = 1 0 0
              imesh  = 4.9 5     iints = 1 1
              jmesh  = 5         jints = 1
              kmesh  = 1         kints = 1
              emesh  = 1e-11 599ilog 10 49i 20 9i 30 14i 60 3i 80
                       7i 160 i 200 19i 1000
              factor = 3.113519e13
c
  fc234       Peak photon spectrum around the chamber (VITAMIN-J-42)
  fmesh234:p  geom   = cyl
              origin = 0 0 53.38
              axs    = 0 0 1
              vec    = 1 0 0
              imesh  = 4.9 5     iints = 1 1
              jmesh  = 5         jints = 1
              kmesh  = 1         kints = 1
              emesh  = 0.001 0.01 i 0.03 i 0.06 0.07 0.075 0.1 i 0.2 i 0.4
                       0.45 0.51 0.512 0.6 i 0.8 1 1.33 1.34 i 1.66
                       2 11i 8 2i 14 20 30 50
              factor = 3.113519e13
c
  fc244       Peak 1-MeV(Si) equivalent flux around the chamber
  fmesh244:n  geom   = cyl
              origin = 0 0 53.38
              axs    = 0 0 1
              vec    = 1 0 0
              imesh  = 4.9 5     iints = 1 1
              jmesh  = 5         jints = 1
              kmesh  = 1         kints = 1
              factor = 3.113519e13
  fm244       10.5263 311 444
c
  fc254       Peak neutron kerma in silicon around the chamber
  fmesh254:n  geom   = cyl
              origin = 0 0 53.38
              axs    = 0 0 1
              vec    = 1 0 0
              imesh  = 4.9 5     iints = 1 1
              jmesh  = 5         jints = 1
              kmesh  = 1         kints = 1
              factor = 3.113519e13
  fm254       3.43541e-10 311 -1 -4
c
  fc264       Peak photon kerma in silicon around the chamber
  fmesh264:p  geom   = cyl
              origin = 0 0 53.38
              axs    = 0 0 1
              vec    = 1 0 0
              imesh  = 4.9 5     iints = 1 1
              jmesh  = 5         jints = 1
              kmesh  = 1         kints = 1
              factor = 3.113519e13
  fm264       3.43541e-10 311 -5 -6
c
c     Neutron flux
c
  fc414       Neutron flux, xy plane at z = 55.88 cm (as plot_xy_room)
  fmesh414:n  geom   = xyz
              origin = -797.56 -482.6 50.8
              imesh  =  624.84  iints = 280
              jmesh  =   482.6  jints = 190
              kmesh  =   60.96  kints =   1
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
  fc424       Neutron flux, xz plane at y = 0 cm (as plot_xz_room)
  fmesh424:n  geom   = xyz
              origin = -797.56 -5.08 -340.36
              imesh  =  624.84  iints = 280
              jmesh  =    5.08  jints =   1
              kmesh  =   177.8  kints = 102
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
  fc434       Neutron flux, yz plane at x = 0 cm (as plot_yz_room)
  fmesh434:n  geom   = xyz
              origin = -5.08 -482.6 -340.36
              imesh  =    5.08  iints =   1
              jmesh  =     381  jints = 170
              kmesh  =   177.8  kints = 102
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
  fc444       Neutron flux, xy plane at z = 55.88 cm (as plot_xy_chamber)
  fmesh444:n  geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
  fc454       Neutron flux, xz plane at y = 0 cm (as plot_xz_chamber)
  fmesh454:n  geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
  fc464       Neutron flux, yz plane at x = 0 cm (as plot_yz_chamber)
  fmesh464:n  geom   = xyz
              origin = -5.08 -76.2 -20.32
              imesh  =    5.08  iints =   1
              jmesh  =    76.2  jints = 120
              kmesh  =  132.08  kints = 120
              emesh  = 5e-7 0.1 10 1000
              factor = 3.113519e13
c
c     Photon flux
c
  fc514       Photon flux, xy plane at z = 55.88 cm (as plot_xy_room)
  fmesh514:p  geom   = xyz
              origin = -797.56 -482.6 50.8
              imesh  =  624.84  iints = 280
              jmesh  =   482.6  jints = 190
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
c
  fc524       Photon flux, xz plane at y = 0 cm (as plot_xz_room)
  fmesh524:p  geom   = xyz
              origin = -797.56 -5.08 -340.36
              imesh  =  624.84  iints = 280
              jmesh  =    5.08  jints =   1
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
c
  fc534       Photon flux, yz plane at x = 0 cm (as plot_yz_room)
  fmesh534:p  geom   = xyz
              origin = -5.08 -482.6 -340.36
              imesh  =    5.08  iints =   1
              jmesh  =     381  jints = 170
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
c
  fc544       Photon flux, xy plane at z = 55.88 cm (as plot_xy_chamber)
  fmesh544:p  geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
c
  fc554       Photon flux, xz plane at y = 0 cm (as plot_xz_chamber)
  fmesh554:p  geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
c
  fc564       Photon flux, yz plane at x = 0 cm (as plot_yz_chamber)
  fmesh564:p  geom   = xyz
              origin = -5.08 -76.2 -20.32
              imesh  =    5.08  iints =   1
              jmesh  =    76.2  jints = 120
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
c
c     1-MeV(Si) equivalent neutron flux
c
  fc614       1-MeV(Si) flux, xy plane at z = 55.88 cm (as plot_xy_room)
  fmesh614:n  geom   = xyz
              origin = -797.56 -482.6 50.8
              imesh  =  624.84  iints = 280
              jmesh  =   482.6  jints = 190
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
  fm614       10.5263 311 444
c
  fc624       1-MeV(Si) flux, xz plane at y = 0 cm (as plot_xz_room)
  fmesh624:n  geom   = xyz
              origin = -797.56 -5.08 -340.36
              imesh  =  624.84  iints = 280
              jmesh  =    5.08  jints =   1
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
  fm624       10.5263 311 444
c
  fc634       1-MeV(Si) flux, yz plane at x = 0 cm (as plot_yz_room)
  fmesh634:n  geom   = xyz
              origin = -5.08 -482.6 -340.36
              imesh  =    5.08  iints =   1
              jmesh  =     381  jints = 170
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
  fm634       10.5263 311 444
c
  fc644       1-MeV(Si) flux, xy plane at z = 55.88 cm (as plot_xy_chamber)
  fmesh644:n  geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
  fm644       10.5263 311 444
c
  fc654       1-MeV(Si) flux, xz plane at y = 0 cm (as plot_xz_chamber)
  fmesh654:n  geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
  fm654       10.5263 311 444
c
  fc664       1-MeV(Si) flux, yz plane at x = 0 cm (as plot_yz_chamber)
  fmesh664:n  geom   = xyz
              origin = -5.08 -76.2 -20.32
              imesh  =    5.08  iints =   1
              jmesh  =    76.2  jints = 120
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
  fm664       10.5263 311 444
c
c     Neutron kerma in silicon [rad(Si)/s]
c
  fc714       Neutron kerma (Si), xy plane at z = 55.88 cm (as plot_xy_room)
  fmesh714:n  geom   = xyz
              origin = -797.56 -482.6 50.8
              imesh  =  624.84  iints = 280
              jmesh  =   482.6  jints = 190
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
  fm714       3.43541e-10 311 -1 -4
c
  fc724       Neutron kerma (Si), xz plane at y = 0 cm (as plot_xz_room)
  fmesh724:n  geom   = xyz
              origin = -797.56 -5.08 -340.36
              imesh  =  624.84  iints = 280
              jmesh  =    5.08  jints =   1
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
  fm724       3.43541e-10 311 -1 -4
c
  fc734       Neutron kerma (Si), yz plane at x = 0 cm (as plot_yz_room)
  fmesh734:n  geom   = xyz
              origin = -5.08 -482.6 -340.36
              imesh  =    5.08  iints =   1
              jmesh  =     381  jints = 170
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
  fm734       3.43541e-10 311 -1 -4
c
  fc744       Neutron kerma (Si), xy plane at z = 55.88 cm (as plot_xy_chamber)
  fmesh744:n  geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
  fm744       3.43541e-10 311 -1 -4
c
  fc754       Neutron kerma (Si), xz plane at y = 0 cm (as plot_xz_chamber)
  fmesh754:n  geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
  fm754       3.43541e-10 311 -1 -4
c
  fc764       Neutron kerma (Si), yz plane at x = 0 cm (as plot_yz_chamber)
  fmesh764:n  geom   = xyz
              origin = -5.08 -76.2 -20.32
              imesh  =    5.08  iints =   1
              jmesh  =    76.2  jints = 120
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
  fm764       3.43541e-10 311 -1 -4
c
c     Photon kerma in silicon [rad(Si)/s]
c
  fc814       Photon kerma (Si), xy plane at z = 55.88 cm (as plot_xy_room)
  fmesh814:p  geom   = xyz
              origin = -797.56 -482.6 50.8
              imesh  =  624.84  iints = 280
              jmesh  =   482.6  jints = 190
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
  fm814       3.43541e-10 311 -5 -6
c
  fc824       Photon kerma (Si), xz plane at y = 0 cm (as plot_xz_room)
  fmesh824:p  geom   = xyz
              origin = -797.56 -5.08 -340.36
              imesh  =  624.84  iints = 280
              jmesh  =    5.08  jints =   1
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
  fm824       3.43541e-10 311 -5 -6
c
  fc834       Photon kerma (Si), yz plane at x = 0 cm (as plot_yz_room)
  fmesh834:p  geom   = xyz
              origin = -5.08 -482.6 -340.36
              imesh  =    5.08  iints =   1
              jmesh  =     381  jints = 170
              kmesh  =   177.8  kints = 102
              factor = 3.113519e13
  fm834       3.43541e-10 311 -5 -6
c
  fc844       Photon kerma (Si), xy plane at z = 55.88 cm (as plot_xy_chamber)
  fmesh844:p  geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              factor = 3.113519e13
  fm844       3.43541e-10 311 -5 -6
c
  fc854       Photon kerma (Si), xz plane at y = 0 cm (as plot_xz_chamber)
  fmesh854:p  geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
  fm854       3.43541e-10 311 -5 -6
c
  fc864       Photon kerma (Si), yz plane at x = 0 cm (as plot_yz_chamber)
  fmesh864:p  geom   = xyz
              origin = -5.08 -76.2 -20.32
              imesh  =    5.08  iints =   1
              jmesh  =    76.2  jints = 120
              kmesh  =  132.08  kints = 120
              factor = 3.113519e13
  fm864       3.43541e-10 311 -5 -6
c
c ------------------------------------------------------------------------------
c     Materials (from the PNNL compendium, ENDF/B-VIII.1)
c ------------------------------------------------------------------------------
c
c         Default data for every material: ENDF/B-VIII.1 neutrons (Lib81, .10c) and
c         photons (eprdata25, .25p), as the OpenMC model's silicon tables
c         (../../data/silicon_*.txt)
c
  m0      nlib=10c  plib=25p
c
c ------------------------------------------------------------------------------
c
c         Air (dry, near sea level)
c         rho = 0.001205 g/cm^3
c
  m4       6012.10c  1.4858312730e-04      6013.10c  1.6070347337e-06
           7014.10c  7.8157536121e-01      7015.10c  2.8553277076e-03
           8016.10c  2.1023586223e-01      8017.10c  8.0084232333e-05
           8018.10c  4.3203335864e-04     18036.10c  1.5718389806e-05
          18038.10c  2.9521611760e-06     18040.10c  4.6524705503e-03
c
c ------------------------------------------------------------------------------
c
c         Concrete, Portland
c         rho = 2.3 g/cm^3
c
  m108     1001.10c  1.6873924210e-01      1002.10c  1.9407244675e-05
           6012.10c  1.4010687221e-03      6013.10c  1.5153578617e-05
           8016.10c  5.6115663647e-01      8017.10c  2.1375895612e-04
           8018.10c  1.1531733159e-03     11023.10c  1.1838220257e-02
          12024.10c  1.1056218375e-03     12025.10c  1.3996984903e-04
          12026.10c  1.5410680379e-04     13027.10c  2.1353800208e-02
          14028.10c  1.8824138727e-01     14029.10c  9.5628086201e-03
          14030.10c  6.3112495738e-03     19039.10c  5.2744001423e-03
          19040.10c  6.6171712339e-07     19041.10c  3.8064004990e-04
          20040.10c  1.8103182585e-02     20042.10c  1.2082358479e-04
          20043.10c  2.5210485234e-05     20044.10c  3.8954868294e-04
          20046.10c  7.4697734026e-07     20048.10c  3.4921190657e-05
          26054.10c  2.4924598409e-04     26056.10c  3.9126289178e-03
          26057.10c  9.0359664721e-05     26058.10c  1.2025212577e-05
c
c ------------------------------------------------------------------------------
c
c         Copper (Cu)
c         rho = 8.96 g/cm^3
c
  m112    29063.10c  6.9150000000e-01     29065.10c  3.0850000000e-01
c
c ------------------------------------------------------------------------------
c
c         Earth, U.S. Average
c         rho = 1.52 g/cm^3
c
  m116     8016.10c  6.6897625085e-01      8017.10c  2.5483021274e-04
           8018.10c  1.3747419371e-03     11023.10c  5.5780533239e-03
          12024.10c  9.0303941357e-03     12025.10c  1.1432325783e-03
          12026.10c  1.2586990687e-03     13027.10c  5.3073574433e-02
          14028.10c  1.8597889604e-01     14029.10c  9.4478723089e-03
          14030.10c  6.2353940617e-03     19039.10c  7.1372817329e-03
          19040.10c  8.9543102717e-07     19041.10c  5.1507947855e-04
          20040.10c  2.5848275793e-02     20042.10c  1.7251559647e-04
          20043.10c  3.5996299110e-05     20044.10c  5.5620948106e-04
          20046.10c  1.0665570107e-06     20048.10c  4.9861540248e-05
          22046.10c  1.6573118377e-04     22047.10c  1.4945939482e-04
          22048.10c  1.4809336809e-03     22049.10c  1.0867947930e-04
          22050.10c  1.0405909478e-04     25055.10c  2.7182020995e-04
          26054.10c  1.2303839620e-03     26056.10c  1.9314396929e-02
          26057.10c  4.4605365534e-04     26058.10c  5.9361553000e-05
c
c ------------------------------------------------------------------------------
c
c         Polyethylene, Non-borated (C2H4)
c         rho = 0.93 g/cm^3
c
c         As HDPE pellets at 75% of 0.95 g/cm^3 (0.7125 g/cm^3) in the
c         beamline penetration, as in the full model
c
  m275     1001.10c  6.6659000000e-01      1002.10c  7.6666666667e-05
           6012.10c  3.2976666667e-01      6013.10c  3.5666666667e-03
  mt275  h-poly.70t
c
c ------------------------------------------------------------------------------
c
c         Silicon (Si)
c         rho = 2.33 g/cm^3
c
c         In no cell: only for its damage energy (MT 444) in the fm cards of the
c         1-MeV(Si) equivalent flux tallies
c
  m311    14028.10c  9.2223000000e-01     14029.10c  4.6850000000e-02
          14030.10c  3.0920000000e-02
c
c ------------------------------------------------------------------------------
c
c         Steel, Stainless 304
c         rho = 8.03 g/cm^3
c
  m331     6012.10c  3.5965407501e-03      6013.10c  3.8899207546e-05
          14028.10c  1.7922298752e-02     14029.10c  9.1046669109e-04
          14030.10c  6.0088858247e-04     15031.10c  7.9296709964e-04
          16032.10c  4.8507331818e-04     16033.10c  3.8299293466e-06
          16034.10c  2.1702932964e-05     16036.10c  5.1065724621e-08
          24050.10c  8.6658055923e-03     24052.10c  1.6711143493e-01
          24053.10c  1.8949095266e-02     24054.10c  4.7168308920e-03
          25055.10c  1.9869795658e-02     26054.10c  3.9042933613e-02
          26056.10c  6.1289056129e-01     26057.10c  1.4154315881e-02
          26058.10c  1.8836796029e-03     28058.10c  6.0141059320e-02
          28060.10c  2.3166228378e-02     28061.10c  1.0070199072e-03
          28062.10c  3.2108201181e-03     28064.10c  8.1770122474e-04
  mt331   fe-56.70t
c
c ------------------------------------------------------------------------------
c
c         Water, Liquid (H2O)
c         rho = 0.997 g/cm^3
c
  m392     1001.10c  6.6659000000e-01      1002.10c  7.6666666667e-05
           8016.10c  3.3252333333e-01      8017.10c  1.2666666667e-04
           8018.10c  6.8333333333e-04
  mt392   h-h2o.70t
c
c ------------------------------------------------------------------------------
