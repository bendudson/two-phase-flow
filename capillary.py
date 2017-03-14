# 
# Plot bubble position over time

from boutdata import collect
from boututils.run_wrapper import shell, launch, getmpirun

from numpy import zeros
import matplotlib.pyplot as plt


def analyse(path="capillary"):
    dx = collect("dx", path=path)[0,0]
    dz = collect("dz", path=path)

    t_array = collect("t_array", path=path)

    vof = collect("vof", path=path)
    nt,nx,_,nz = vof.shape

    mid_x = (nx-4.)/2. + 1.5
    mid_z = nz / 2.0 - 0.5  # Index at the middle

    imz = int(nz/4)

    result = zeros(nt)

    # Find where vof goes from 0 to 1 or 1 to 0
    for tind in range(nt):
        data = vof[tind, :, 0,imz]
    
        indl = 0
        indh = nx-1

        if data[indl] < 0.5 and data[indh] > 0.5:
            # Find where vof transitions from 0 to 1
            data = 1.0 - data
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
    return t_array, (mid_x - result) * dx # meters

MPIRUN=getmpirun()

##########################################################
# lambda = 0.1m test case

cases = [
    (r"FD, $D=0.1$", "model:curv_method=0 model:vof_D=0.1", '-ko')
    ,(r"HF, $D=0.1$", "model:curv_method=1 model:vof_D=0.1", '-rx')
    ,(r"FD, $D=1$", "model:curv_method=0 model:vof_D=1", '--ko')
    ,(r"HF, $D=1$", "model:curv_method=1 model:vof_D=1", '--rx')
]

for label, opts, sym in cases:
    print("Running case: "+label)
    s, out = launch("./two-phase-flow -d capillary "+opts, runcmd=MPIRUN, nproc=1, pipe=True)
    time, disp = analyse("capillary")
    plt.plot(time, disp, sym, label=label)

plt.grid()
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.legend()
plt.xlim([0,0.6])
plt.ylim([-0.01, 0.01])
plt.title(r"Capillary wave $\lambda = 0.1$m")
plt.savefig("capillary_1.pdf")
plt.savefig("capillary_1.png")

plt.show()

##########################################################
# lambda = 1m test case

for label, opts, sym in cases:
    print("Running case: "+label)
    s, out = launch("./two-phase-flow -d capillary mesh:Lz=1 nout=100 timestep=0.2 "+opts, runcmd=MPIRUN, nproc=1, pipe=True)
    time, disp = analyse("capillary")
    plt.plot(time, disp, sym, label=label)

plt.grid()
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.legend()
plt.xlim([0,20])
plt.ylim([-0.1, 0.1])
plt.title(r"Capillary wave $\lambda = 1$m")
plt.savefig("capillary_2.pdf")
plt.savefig("capillary_2.png")

plt.show()
