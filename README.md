Stokes (up to) 5th-order Regular Wave
Implemented according to Fenton, 1985

LM, 16 April 2015
[git SHA1 ID: 3e2e3bb61a7dcb3769801399c4bfc66f99a7620d]

Functions:
fDispersionV5.m - dispersion equation, variable orders from 1-5
fStokesIn.m - (optional) use this prepare your inputs
fStokesEta.m - surface profile
fStokesVel.m - velocities
fStokesAcc.m - accelerations
fFitStokes.m - interactive tool to fit a stokes wave to an exisiting wave profile

Input preparation:
1) to use each function individually, read function helps
2) for the input wave properties, all functions (except for fDispersionV5.m) require a wave properties struct (let's call it 'wp'). The format of this struct is:
   wp.H - wave height
   wp.T - wave period
   wp.d - depth from bed to SWL
   wp.order - from 1 up to 5
   wp.omega - circular wave frequency
   wp.k - wave number, may be calculated with fDispersionV5.m
   wp.lamda - wave length
   wp.ReturnFlow: [true | false] - switch to incorporated Eularien return flow, will affect your wave number, and horizontal velocity. Always OFF for 1st-order calculations.
   wp.SwlAdjust: [true | false] - switch to shift the SWL (according to Bowden), correct up to 3rd-order. Only affects surface profile. Always OFF for 1st-order calculations.
   wp.DTerms: [true | false] - switch to include D coefficients in Fenton 1985 in the calculation of wave number k. Based on the effect that the wave is propagating on top of stokes drift. Will affect your wave number.
3) you are encouraged to use the fStokesIn.m function to prepare the above wp struct for you. If you choose to do this yourself please ensure all settings are consistent throughout.

Defaults - ignore these if you always specify all options:
1) fStokesEta.m, fStokesVel.m & fStokesAcc.m expect all fields of the wp struct to be present and will report error if not. No default is set for these 3 functions.
2) ReturnFlow, SwlAdjust, DTerms all default to be OFF (if not given), when fStokesIn.m is used to prepare the inputs. Same goes for fDispersionV5.m.
3) fDispersionV5.m defaults to 1st order when H is not given; and defaults to 5th order when H is given.

