#!/usr/bin/env python
#
# Reads data, produces sequence of images for animation

from boutdata import collect
from numpy import transpose, amax, amin, linspace, arange
import matplotlib.pyplot as plt


path="data"  # The path to the data to analyse


dx = collect("dx", path=path)[0,0]
dz = collect("dz", path=path)

vort = collect("vorticity",path=path)
vof = collect("vof",path=path)

nt,nx,_,nz = vort.shape

x = (arange(nx)-.5)*dx
z = (arange(nz)+0.5)*dz

#vort = transpose(vort[:,:,0,:], axes=[0,2,1])
#vof = transpose(vof[:,:,0,:], axes=[0,2,1])
vort = vort[:,:,0,:]
vof = vof[:,:,0,:]

def calc_levels(var, num=50):
  ma = amax(var)
  mi = amin(var)
  return linspace(mi, ma, num=num)

vort_lev = calc_levels(vort)
vof_lev = calc_levels(vof)

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

for tind in range(nt):
  ax1.clear()
  ax2.clear()

  ax1.contourf(z,x,vort[tind,:,:], vort_lev)
  ax1.set_title("Vorticity")
  ax1.set_xlim([0,1])
  ax1.set_ylim([0,2])

  ax2.contourf(z,x,vof[tind,:,:], vof_lev)
  ax2.set_title("Fluid type")
  ax2.set_xlim([0,1])
  ax2.set_ylim([0,2])

  plt.savefig("two-phase-%04d.png" % (tind))
plt.show()

