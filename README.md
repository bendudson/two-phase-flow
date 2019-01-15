Incompressible three phase flow in 2D
====================================

A simple model for three fluids (e.g. liquid, gas, solid) in 2D,
with different densities and viscosities. Includes models
for buoyancy and surface tension.

The flow is assumed to be incompressible, so a vorticity equation
is evolved for the momentum equation and no pressure equation is
needed. Uses Arakawa's 2nd-order method to advect the vorticity,
and upwinding with regularised anti-diffusion to advect the fluid type
and track the location of boundaries between fluids.

Compiled with [BOUT++ version 4.0](https://github.com/boutproject/BOUT-dev/releases/tag/v4.0.0)

Stationary bubble
-----------------

A test of the handling of surface tension is a stationary bubble.
This should produce no flows, but numerically parasitic flows are observed.

$ ./two-phase-flow -d bubble
$ python bubble.py


Capillary wave
--------------

A boundary between two viscous fluids is initially perturbed. The
frequency and damping rate of the resulting oscillation can be compared
to analytic theory and previous results

$ python capillary.py

