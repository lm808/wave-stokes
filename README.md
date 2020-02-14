# wave-stokes

An implimentation of the Stokes 5th-order regular wave theory, as proposed by Fenton (1985).

Released under the GNU General Public License v3.0. The author reserves the rights to later publications of the project under a different licence. All usage must incorporate appropriate citation as listed below.

Copyright (c) 2013-2020 lm808. All rights reserved.

# Cite as
L. Ma, wave-stokes (2020), open-source repository, https://github.com/lm808/wave-stokes
Fenton, J. (1985). A 5 th -order stokes theory for steady waves. Journal of Waterway, Port, Coastal, and Ocean Engineering, 111 , 216–234.

# Usage
1. Download and add to MATLAB path.
2. Type "help <function name>" to see further documentation.

# List of functions
`fStokesIn.m`
Function to prepare the input. This will output a MATLAB containing all the necessary inputs. Read the comments inside the function for details.

`fDispersionV5.m`
The dispersion equation linking the wave period and wave length, with variable orders (in terms of the wave steepness, Hk/2) from 1-5. If the order is not given, it defaults to 1st-order when H is not given; and defaults to 5th-order when H is given.

`fStokesEta.m`
Computes the free water surface elevation.

`fStokesVel.m`
Computes the underlying water particle velocities.

`fStokesAcc.m`
Computes the underlying water particle accelerations.

`fStokesFit.m`
An interactive tool to fit a stokes wave to an exisiting wave profile.
