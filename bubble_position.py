# 
# Plot bubble position over time

from boutdata import collect

from numpy import zeros
import matplotlib.pyplot as plt

path="bubble/"

dx = collect("dx", path=path)
dz = collect("dz", path=path)

vof = collect("vof", path=path)
nt,nx,_,nz = vof.shape

mid_x = (nx-4.)/2. + 1.5
mid_z = nz / 2.0 - 0.5  # Index at the middle

imx = int(mid_x)
imz = int(mid_z)

result = zeros(nt)

for tind in range(nt):
    # Find where vof transitions from 1 to 0
    data = vof[tind, :, 0,imz]
    
    indl = imx
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
    
plt.plot(result)
plt.show()
