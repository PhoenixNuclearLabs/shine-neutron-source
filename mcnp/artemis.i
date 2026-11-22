SHINE Artemis neutron generator, simplified model
c
c  Author:       Lucas Jacobson
c  Created on:   November 20, 2026
c  Last updated: November 21, 2026
c
c  Copyright 2026, SHINE Technologies. All rights reserved.
c
c  A very simplified model of the Artemis neutron generator below ground:
c  the target chamber with its double wall, the upstream section of the
c  beamline (the beam tube runs up to the ground), the two irradiation
c  cavities with their cassette regions and ducts, the water pit (modeled
c  full of water up to the ground), its liner, the concrete around it, and
c  the earth. The cavity and duct walls are all 1.27 cm thick. The model
c  is symmetric about x = 0 and y = 0, and every cell is a plain
c  intersection of half-spaces.
c
c  The backstop is at the origin and the beam travels down (-z) from the
c  top of the 137 cm source region. The ground is at z = 289.56 cm.
c  Tallies are scaled by the total neutron yield (yield_total in
c  ../source/Artemis/summary.toml).
c
c  Geometry plots: mcnp6 ip i=artemis.i com=plot_xy_pit (or plot_xz_pit, or
c  the zoomed plot_xy_cavity and plot_xz_cavity). The mesh tallies cover the
c  same windows as these plots.
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
    1    0                    -103   412  -420                    $ Vacuum in target chamber
    2  112  -8.96             -103   411  -412                    $ Copper beam stop
    3  331  -8.03             -104   410  -411                    $ Steel plate below beam stop
    4  331  -8.03        103  -104   411  -419                    $ Target chamber inner wall
    5  392  -0.997       104  -105   410  -419                    $ Cooling water jacket
    6  392  -0.997            -105   409  -410                    $ Cooling water behind beam stop
    7  331  -8.03        105  -106   409  -419                    $ Target chamber outer wall
    8  331  -8.03             -106   407  -409                    $ Target chamber end plate
    9  331  -8.03        103  -110   419  -420                    $ Upstream tube bottom plate
   10    0                    -109   420  -421                    $ Vacuum in upstream tube
   11  331  -8.03        109  -110   420  -421                    $ Upstream tube
   12  331  -8.03        102  -110   421  -422                    $ Upstream tube top plate
   13    0                    -101   421  -424                    $ Vacuum in beam tube
   14  331  -8.03        101  -102   421  -424                    $ Beam tube
   15  392  -0.997       106  -107   407  -418                    $ Water between chamber and cavities
   16  392  -0.997            -107   406  -407                    $ Water below chamber, in cavities
c
c ------------------------------------------------------------------------------
c     Irradiation cavity, cassette region, and ducts (-x)
c ------------------------------------------------------------------------------
c
  101    4  -0.001205    211  -212   307  -314   408  -417        $ Air in cavity, beyond cassette (-x)
  102    4  -0.001205    212  -213   311  -314   408  -417        $ Air in cavity, +y of cassette (-x)
  103    4  -0.001205    212  -213   307  -310   408  -417        $ Air in cavity, -y of cassette (-x)
  104    4  -0.001205    212  -213   310  -311   408  -415   108  $ Air in cavity, below cassette (-x)
                         vol=38386.77
  105    4  -0.001205    212  -213   310  -311   416  -417   108  $ Air in cavity, above cassette (-x)
                         vol=15344.71
  106    4  -0.001205    212  -213   310  -311   415  -416   108  $ Cassette region (-x)
                         vol=153697.03
  107   13  -2.7         210  -214   306  -315   406  -408   107  $ Cavity bottom (-x)
                         vol=7896.08
  108   13  -2.7         210  -214   306  -315   417  -418   107  $ Cavity top (-x)
                         vol=7896.08
  109   13  -2.7         210  -214   314  -315   408  -417        $ Cavity +y wall (-x)
  110   13  -2.7         210  -214   306  -307   408  -417        $ Cavity -y wall (-x)
  111   13  -2.7         213  -214   307  -314   408  -417   107  $ Cavity inner wall (-x)
                         vol=10332.27
  112   13  -2.7        -213   107  -108   408  -417              $ Cavity wall around chamber (-x)
                         vol=1285.58
  113   13  -2.7         210  -211   312  -314   408  -417        $ Cavity outer wall, +y of duct (-x)
  114   13  -2.7         210  -211   307  -309   408  -417        $ Cavity outer wall, -y of duct (-x)
  115   13  -2.7         210  -211   309  -312   408  -414        $ Cavity outer wall, below duct (-x)
  116    4  -0.001205    208  -211   309  -312   414  -417        $ Air in horizontal duct (-x)
  117   13  -2.7         209  -210   312  -313   413  -418        $ Horizontal duct +y wall (-x)
  118   13  -2.7         209  -210   308  -309   413  -418        $ Horizontal duct -y wall (-x)
  119   13  -2.7         209  -210   309  -312   413  -414        $ Horizontal duct bottom (-x)
  120   13  -2.7         209  -210   309  -312   417  -418        $ Horizontal duct top (-x)
  121    4  -0.001205    207  -208   309  -312   408  -424        $ Air in vertical duct (-x)
  122   13  -2.7         208  -209   309  -312   408  -414        $ Vertical duct inner wall, below (-x)
  123   13  -2.7         208  -209   309  -312   417  -424        $ Vertical duct inner wall, above (-x)
  124   13  -2.7         206  -207   309  -312   408  -424        $ Vertical duct outer wall (-x)
  125   13  -2.7         206  -209   312  -313   406  -424        $ Vertical duct +y wall (-x)
  126   13  -2.7         206  -209   308  -309   406  -424        $ Vertical duct -y wall (-x)
  127   13  -2.7         206  -209   309  -312   406  -408        $ Vertical duct bottom (-x)
c
c ------------------------------------------------------------------------------
c     Irradiation cavity, cassette region, and ducts (+x)
c ------------------------------------------------------------------------------
c
  201    4  -0.001205    216  -217   307  -314   408  -417        $ Air in cavity, beyond cassette (+x)
  202    4  -0.001205    215  -216   311  -314   408  -417        $ Air in cavity, +y of cassette (+x)
  203    4  -0.001205    215  -216   307  -310   408  -417        $ Air in cavity, -y of cassette (+x)
  204    4  -0.001205    215  -216   310  -311   408  -415   108  $ Air in cavity, below cassette (+x)
                         vol=38386.77
  205    4  -0.001205    215  -216   310  -311   416  -417   108  $ Air in cavity, above cassette (+x)
                         vol=15344.71
  206    4  -0.001205    215  -216   310  -311   415  -416   108  $ Cassette region (+x)
                         vol=153697.03
  207   13  -2.7         214  -218   306  -315   406  -408   107  $ Cavity bottom (+x)
                         vol=7896.08
  208   13  -2.7         214  -218   306  -315   417  -418   107  $ Cavity top (+x)
                         vol=7896.08
  209   13  -2.7         214  -218   314  -315   408  -417        $ Cavity +y wall (+x)
  210   13  -2.7         214  -218   306  -307   408  -417        $ Cavity -y wall (+x)
  211   13  -2.7         214  -215   307  -314   408  -417   107  $ Cavity inner wall (+x)
                         vol=10332.27
  212   13  -2.7         215   107  -108   408  -417              $ Cavity wall around chamber (+x)
                         vol=1285.58
  213   13  -2.7         217  -218   312  -314   408  -417        $ Cavity outer wall, +y of duct (+x)
  214   13  -2.7         217  -218   307  -309   408  -417        $ Cavity outer wall, -y of duct (+x)
  215   13  -2.7         217  -218   309  -312   408  -414        $ Cavity outer wall, below duct (+x)
  216    4  -0.001205    217  -220   309  -312   414  -417        $ Air in horizontal duct (+x)
  217   13  -2.7         218  -219   312  -313   413  -418        $ Horizontal duct +y wall (+x)
  218   13  -2.7         218  -219   308  -309   413  -418        $ Horizontal duct -y wall (+x)
  219   13  -2.7         218  -219   309  -312   413  -414        $ Horizontal duct bottom (+x)
  220   13  -2.7         218  -219   309  -312   417  -418        $ Horizontal duct top (+x)
  221    4  -0.001205    220  -221   309  -312   408  -424        $ Air in vertical duct (+x)
  222   13  -2.7         219  -220   309  -312   408  -414        $ Vertical duct inner wall, below (+x)
  223   13  -2.7         219  -220   309  -312   417  -424        $ Vertical duct inner wall, above (+x)
  224   13  -2.7         221  -222   309  -312   408  -424        $ Vertical duct outer wall (+x)
  225   13  -2.7         219  -222   312  -313   406  -424        $ Vertical duct +y wall (+x)
  226   13  -2.7         219  -222   308  -309   406  -424        $ Vertical duct -y wall (+x)
  227   13  -2.7         219  -222   309  -312   406  -408        $ Vertical duct bottom (+x)
c
c ------------------------------------------------------------------------------
c     Water in pit
c ------------------------------------------------------------------------------
c
  301  392  -0.997       205  -223   305  -316   405  -406        $ Water below everything
  302  392  -0.997       205  -223   305  -306   406  -424        $ Water, -y side
  303  392  -0.997       205  -223   315  -316   406  -424        $ Water, +y side
  304  392  -0.997       205  -206   306  -315   406  -424        $ Water, -x end
  305  392  -0.997       222  -223   306  -315   406  -424        $ Water, +x end
  306  392  -0.997       206  -209   306  -308   406  -424        $ Water, -y of vertical duct (-x)
  307  392  -0.997       206  -209   313  -315   406  -424        $ Water, +y of vertical duct (-x)
  308  392  -0.997       219  -222   306  -308   406  -424        $ Water, -y of vertical duct (+x)
  309  392  -0.997       219  -222   313  -315   406  -424        $ Water, +y of vertical duct (+x)
  310  392  -0.997       209  -210   306  -308   406  -418        $ Water, -y of horizontal duct (-x)
  311  392  -0.997       209  -210   313  -315   406  -418        $ Water, +y of horizontal duct (-x)
  312  392  -0.997       209  -210   308  -313   406  -413        $ Water below horizontal duct (-x)
  313  392  -0.997       218  -219   306  -308   406  -418        $ Water, -y of horizontal duct (+x)
  314  392  -0.997       218  -219   313  -315   406  -418        $ Water, +y of horizontal duct (+x)
  315  392  -0.997       218  -219   308  -313   406  -413        $ Water below horizontal duct (+x)
  316  392  -0.997       209  -219   306  -315   418  -419   106  $ Water around target chamber
                         vol=604207.22
  317  392  -0.997       209  -219   306  -315   419  -422   110  $ Water around upstream tube
                         vol=1817972.44
  318  392  -0.997       209  -219   306  -315   422  -424   102  $ Water around beam tube
                         vol=1894446.11
c
c ------------------------------------------------------------------------------
c     Pit liner and air gap
c ------------------------------------------------------------------------------
c
  401  331  -8.03        204  -224   304  -317   404  -405        $ Pit liner, bottom
  402  331  -8.03        204  -224   304  -305   405  -424        $ Pit liner, -y side
  403  331  -8.03        204  -224   316  -317   405  -424        $ Pit liner, +y side
  404  331  -8.03        204  -205   305  -316   405  -424        $ Pit liner, -x side
  405  331  -8.03        223  -224   305  -316   405  -424        $ Pit liner, +x side
  411    4  -0.001205    203  -225   303  -318   403  -404        $ Air gap, bottom
  412    4  -0.001205    203  -225   303  -304   404  -424        $ Air gap, -y side
  413    4  -0.001205    203  -225   317  -318   404  -424        $ Air gap, +y side
  414    4  -0.001205    203  -204   304  -317   404  -424        $ Air gap, -x side
  415    4  -0.001205    224  -225   304  -317   404  -424        $ Air gap, +x side
c
c ------------------------------------------------------------------------------
c     Concrete
c ------------------------------------------------------------------------------
c
  501  108  -2.3         202  -226   302  -319   402  -403        $ Slab below pit
  502  108  -2.3         202  -226   302  -303   403  -423        $ Pit wall, -y side
  503  108  -2.3         202  -226   318  -319   403  -423        $ Pit wall, +y side
  504  108  -2.3         202  -203   303  -318   403  -423        $ Pit wall, -x side
  505  108  -2.3         225  -226   303  -318   403  -423        $ Pit wall, +x side
  506  108  -2.3         201  -227   301  -303   423  -424        $ Floor, -y side
  507  108  -2.3         201  -227   318  -320   423  -424        $ Floor, +y side
  508  108  -2.3         201  -203   303  -318   423  -424        $ Floor, -x side
  509  108  -2.3         225  -227   303  -318   423  -424        $ Floor, +x side
c
c ------------------------------------------------------------------------------
c     Earth
c ------------------------------------------------------------------------------
c
  601  116  -1.52        201  -227   301  -320   401  -402        $ Earth below slab
  602  116  -1.52        201  -227   301  -302   402  -423        $ Earth beside pit walls, -y side
  603  116  -1.52        201  -227   319  -320   402  -423        $ Earth beside pit walls, +y side
  604  116  -1.52        201  -202   302  -319   402  -423        $ Earth beside pit walls, -x side
  605  116  -1.52        226  -227   302  -319   402  -423        $ Earth beside pit walls, +x side
c
c ------------------------------------------------------------------------------
c     Air above ground
c ------------------------------------------------------------------------------
c
  701    4  -0.001205    201  -227   301  -320   424  -425        $ Air above ground
c
c ------------------------------------------------------------------------------
c     Outside the model
c ------------------------------------------------------------------------------
c
  999    0             (-201:227:-301:320:-401:425)               $ Rest of universe

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
  102  cz       4.59740  $ Beam tube
  103  cz       5.44830  $ Target chamber inner radius
  104  cz       5.72770  $ Target chamber inner wall
  105  cz       5.86740  $ Cooling water jacket
  106  cz       6.09600  $ Target chamber outer wall
  107  cz       6.40969  $ Water around the target chamber in the cavities
  108  cz       7.06501  $ Cavity wall around the target chamber
  109  cz      12.06500  $ Upstream tube
  110  cz      12.70000  $ Upstream tube
c
c ------------------------------------------------------------------------------
c     Planes normal to x
c ------------------------------------------------------------------------------
c
  201  px    -347.98000  $ Floor and earth
  202  px    -256.54000  $ Pit walls and slab
  203  px    -215.90000  $ Air gap
  204  px    -195.58000  $ Pit liner
  205  px    -195.23710  $ Pit
  206  px    -190.50000  $ Vertical duct outer wall
  207  px    -189.23000  $ Vertical duct outer wall
  208  px    -128.27000  $ Vertical duct inner wall
  209  px    -127.00000  $ Vertical duct inner wall
  210  px     -67.94500  $ Cavity outer wall
  211  px     -66.67500  $ Cavity outer wall
  212  px     -51.18100  $ Cassette region
  213  px      -1.27000  $ Cavity inner wall
  214  px       0.00000  $ Center
  215  px       1.27000  $ Cavity inner wall
  216  px      51.18100  $ Cassette region
  217  px      66.67500  $ Cavity outer wall
  218  px      67.94500  $ Cavity outer wall
  219  px     127.00000  $ Vertical duct inner wall
  220  px     128.27000  $ Vertical duct inner wall
  221  px     189.23000  $ Vertical duct outer wall
  222  px     190.50000  $ Vertical duct outer wall
  223  px     195.23710  $ Pit
  224  px     195.58000  $ Pit liner
  225  px     215.90000  $ Air gap
  226  px     256.54000  $ Pit walls and slab
  227  px     347.98000  $ Floor and earth
c
c ------------------------------------------------------------------------------
c     Planes normal to y
c ------------------------------------------------------------------------------
c
  301  py    -347.98000  $ Floor and earth
  302  py    -256.54000  $ Pit walls and slab
  303  py    -215.90000  $ Air gap
  304  py    -195.58000  $ Pit liner
  305  py    -195.23710  $ Pit
  306  py     -46.22800  $ Cavities (outside)
  307  py     -44.95800  $ Cavities (inside)
  308  py     -23.81250  $ Ducts (outside)
  309  py     -22.54250  $ Ducts (inside)
  310  py     -20.32000  $ Cassette region
  311  py      20.32000  $ Cassette region
  312  py      22.54250  $ Ducts (inside)
  313  py      23.81250  $ Ducts (outside)
  314  py      44.95800  $ Cavities (inside)
  315  py      46.22800  $ Cavities (outside)
  316  py     195.23710  $ Pit
  317  py     195.58000  $ Pit liner
  318  py     215.90000  $ Air gap
  319  py     256.54000  $ Pit walls and slab
  320  py     347.98000  $ Floor and earth
c
c ------------------------------------------------------------------------------
c     Planes normal to z
c ------------------------------------------------------------------------------
c
  401  pz    -299.72000  $ Bottom of earth
  402  pz    -208.28000  $ Bottom of slab
  403  pz    -177.80000  $ Top of slab
  404  pz    -167.64000  $ Bottom of pit liner
  405  pz    -167.29710  $ Bottom of pit
  406  pz      -4.26720  $ Bottom of cavities and vertical ducts
  407  pz      -3.52552  $ Bottom of target chamber
  408  pz      -2.99720  $ Bottom of air in cavities and vertical ducts
  409  pz      -1.18872  $ Top of end plate (bottom of cooling water)
  410  pz      -0.95250  $ Bottom of steel plate below beam stop
  411  pz      -0.63500  $ Bottom of beam stop
  412  pz       0.00000  $ Backstop (top of beam stop)
  413  pz       9.70280  $ Bottom of horizontal ducts
  414  pz      10.97280  $ Bottom of air in horizontal ducts
  415  pz      16.51000  $ Bottom of cassette regions
  416  pz      94.61500  $ Top of cassette regions
  417  pz     102.41280  $ Top of air in cavities and horizontal ducts
  418  pz     103.68280  $ Top of cavities and horizontal ducts
  419  pz     129.54000  $ Top of target chamber walls
  420  pz     130.17500  $ Top of upstream tube bottom plate
  421  pz     208.02600  $ Top of upstream tube
  422  pz     208.66100  $ Top of upstream tube top plate
  423  pz     269.24000  $ Bottom of floor
  424  pz     289.56000  $ Ground
  425  pz     751.84000  $ Top of air above ground

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
  nps         1e8
  prdmp       j 1e7 1 2 j
c
  mode        n p
  imp:n       1 112r 0
  imp:p       1 112r 0
c
c ------------------------------------------------------------------------------
c     Source definition
c ------------------------------------------------------------------------------
c
  read        file = ../source/Artemis/source.sdef
c
c ------------------------------------------------------------------------------
c     Tallies
c ------------------------------------------------------------------------------
c
  fc14        Neutron energy spectrum in the cassette regions (CCFE-709)
  f14:n       106 206
  e14         1e-11 599ilog 10 49i 20 9i 30 14i 60 3i 80 7i 160 i 200 19i 1000
  fq14        e f
  fm14        1.77562e13
c
  fc24        Photon energy spectrum in the cassette regions (VITAMIN-J-42)
  f24:p       106 206
  e24         0.001 0.01 i 0.03 i 0.06 0.07 0.075 0.1 i 0.2 i 0.4 0.45 0.51
              0.512 0.6 i 0.8 1 1.33 1.34 i 1.66 2 11i 8 2i 14 20 30 50
  fq24        e f
  fm24        1.77562e13
c
c ------------------------------------------------------------------------------
c     Mesh tallies
c ------------------------------------------------------------------------------
c
c     Neutrons
c
  fc34        Neutron flux in the xy plane at z = 55.88 cm (as plot_xy_pit)
  fmesh34:n   geom   = xyz
              origin = -243.84 -243.84 50.8
              imesh  =  243.84  iints = 192
              jmesh  =  243.84  jints = 192
              kmesh  =   60.96  kints =   1
              factor = 1.77562e13
c
  fc44        Neutron flux in the xz plane at y = 0 cm (as plot_xz_pit)
  fmesh44:n   geom   = xyz
              origin = -243.84 -5.08 -187.96
              imesh  =  243.84  iints = 192
              jmesh  =    5.08  jints =   1
              kmesh  =  299.72  kints = 192
              factor = 1.77562e13
c
  fc54        Neutron flux in the xy plane at z = 55.88 cm (as plot_xy_cavity)
  fmesh54:n   geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              factor = 1.77562e13
c
  fc64        Neutron flux in the xz plane at y = 0 cm (as plot_xz_cavity)
  fmesh64:n   geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              factor = 1.77562e13
c
c     Photons
c
  fc74        Photon flux in the xy plane at z = 55.88 cm (as plot_xy_pit)
  fmesh74:p   geom   = xyz
              origin = -243.84 -243.84 50.8
              imesh  =  243.84  iints = 192
              jmesh  =  243.84  jints = 192
              kmesh  =   60.96  kints =   1
              factor = 1.77562e13
c
  fc84        Photon flux in the xz plane at y = 0 cm (as plot_xz_pit)
  fmesh84:p   geom   = xyz
              origin = -243.84 -5.08 -187.96
              imesh  =  243.84  iints = 192
              jmesh  =    5.08  jints =   1
              kmesh  =  299.72  kints = 192
              factor = 1.77562e13
c
  fc94        Photon flux in the xy plane at z = 55.88 cm (as plot_xy_cavity)
  fmesh94:p   geom   = xyz
              origin = -76.2 -76.2 50.8
              imesh  =    76.2  iints = 120
              jmesh  =    76.2  jints = 120
              kmesh  =   60.96  kints =   1
              factor = 1.77562e13
c
  fc104       Photon flux in the xz plane at y = 0 cm (as plot_xz_cavity)
  fmesh104:p  geom   = xyz
              origin = -76.2 -5.08 -20.32
              imesh  =    76.2  iints = 120
              jmesh  =    5.08  jints =   1
              kmesh  =  132.08  kints = 120
              factor = 1.77562e13
c
c ------------------------------------------------------------------------------
c     Materials (from the PNNL compendium, ENDF/B-VIII.1)
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
c         Aluminum, alloy 6061-O
c         rho = 2.7 g/cm^3
c
  m13     12024.10c  8.8169277578e-03     12025.10c  1.1162080969e-03
          12026.10c  1.2289451147e-03     13027.10c  9.7733008842e-01
          14028.10c  5.3450290793e-03     14029.10c  2.7153162700e-04
          14030.10c  1.7920507805e-04     22046.10c  4.0960555374e-05
          22047.10c  3.6938973573e-05     22048.10c  3.6601359299e-04
          22049.10c  2.6860194494e-05     22050.10c  2.5718263859e-05
          24050.10c  4.4207409019e-05     24052.10c  8.5249587900e-04
          24053.10c  9.6666189433e-05     24054.10c  2.4062260605e-05
          25055.10c  4.3258616015e-04     26054.10c  1.1607847426e-04
          26056.10c  1.8221838027e-03     26057.10c  4.2082170564e-05
          26058.10c  5.6003643696e-06     29063.10c  8.1185431728e-04
          29065.10c  3.6219386389e-04     30064.10c  2.9229283659e-04
          30066.10c  1.6940606828e-04     30067.10c  2.4840167724e-05
          30068.10c  1.1520218205e-04     30070.10c  3.8210984480e-06
  mt13    al-27.70t
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
