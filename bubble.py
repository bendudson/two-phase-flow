#!/usr/bin/env python
# 

path="bubble"

from boutdata import collect
from numpy import arange, mean, sqrt, roll

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size':14})

vort = collect("vort", path=path)
dx = collect("dx",path=path)[0,0]
dz = collect("dz", path=path)

time = collect("t_array", path=path) * 1e3

nt,nx,_,nz = vort.shape

mid_x = (nx-4.)/2. + 1.5

x = (arange(nx)-(nx/2.)-0.5) * dx * 1e3
z = (arange(nz)-(nz/2.)+0.5) * dz * 1e3

# Calculate RMS vorticity
vrms = sqrt(mean(mean(vort**2,axis=3),axis=1))

plt.plot(time, vrms)
plt.xlim([0,100])
plt.xlabel("Time [ms]")
plt.ylabel(r"RMS vorticity [kg/m$^3$/s]")
plt.savefig("bubble_vorticity.pdf")
plt.savefig("bubble_vorticity.png")
plt.show()

vof = collect("vof", path=path, tind=nt-1)
psi = collect("psi", path=path, tind=nt-1)

plt.contourf(z, x, vof[-1,:,0,:], 50)
plt.contour(z,x,psi[-1,:,0,:], 10)
plt.xlabel("X [mm]")
plt.ylabel("Z [mm]")

plt.savefig("bubble_streamlines.pdf")
plt.savefig("bubble_streamlines.png")

plt.show()

# Calculate velocities
psi = psi[-1,:,0,:]

vx = (roll(psi, 1, axis=1) - roll(psi, -1, axis=1))/(2.*dz)
vz = -(roll(psi, 1, axis=0) - roll(psi, -1, axis=0))/(2.*dx)

vmag = sqrt(vx**2 + vz**2)

plt.contourf(z,x,vmag,50)
plt.colorbar()
plt.quiver(z,x,vz,vx, scale=3e-2)
plt.contour(z,x,vof[-1,:,0,:], 5, colors='k')
#plt.contour(z,x,psi, 20)
plt.xlabel("X [mm]")
plt.ylabel("Z [mm]")
plt.title("Flows around a stationary bubble [m/s]")

plt.savefig("bubble_flows.pdf")
plt.savefig("bubble_flows.png")

plt.show()
