# 
# Plot bubble position over time

from boutdata import collect

from numpy import zeros
import matplotlib.pyplot as plt

path="capillary/"

dx = collect("dx", path=path)[0,0]
dz = collect("dz", path=path)

vof = collect("vof", path=path)
nt,nx,_,nz = vof.shape

mid_x = (nx-4.)/2. + 1.5
mid_z = nz / 2.0 - 0.5  # Index at the middle

imz = int(nz/4)

result = zeros(nt)

for tind in range(nt):
    # Find where vof transitions from 1 to 0
    data = 1.0 - vof[tind, :, 0,imz]
    
    indl = 0
    indh = nx-1
    while indl != indh:
        mid = int((indl + indh)/2.0)
        
        if (mid == indl) or (mid == indh):
            break

        if data[mid] > 0.5:
            indl = mid
        else:
            indh = mid
    
    if data[mid] < 0.5:
        mid -= 1
    mid += (0.5 - data[mid]) / (data[mid+1] - data[mid])
    
    result[tind] = mid
    
# convert to displacement
result = (mid_x - result) * dx # meters

print nx*dx
print nz*dz

plt.plot(result)
plt.show()
