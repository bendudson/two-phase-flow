Incompressible two-phase flow in 2D
===================================

A simple model for two fluids (e.g. liquid and gas) in 2D,
with different densities and viscosities. Buoyancy is included,
but not surface tension (yet). 

The flow is assumed to be incompressible, so a vorticity equation
is evolved for the momentum equation and no pressure equation is
needed. Uses Arakawa's 2nd-order method to advect the vorticity,
and WENO3 for the fluid type field.

Compiled with [BOUT++ version 4.0](https://github.com/boutproject/BOUT-dev/releases/tag/v4.0.0)

