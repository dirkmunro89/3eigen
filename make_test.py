import json
import numpy as np

out = {}

ex = 50;
ey = 50;
ez = 50;
e = ex*ey*ez
gss = 2;
mat = [1.0, .3, "plane_stress"]

lx=1.
ly=1.
lz=1.

nx = ex+1;
ny = ey+1;
nz = ez+1;
n = nx*ny*nz;

nds = []
els = []
bds = []
lds = []

n=0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            nds.append([float(i)*lx/ex,float(j)*ly/ey,float(k)*lz/ez] )
            for d in range(3):
                if k == 0:
                    bds.append([n,0.])
                if k == nz-1 and j == ny-1 and d == 1:
                    if i == 0 or i == nx-1:
                        lds.append([n,1.])
                    else:
                        lds.append([n,1.])
                n=n+1

n1=0
n2=1
n3=nx+1
n4=nx
n5=nx*ny
n6=nx*ny+1
n7=nx*ny+nx+1
n8=nx*ny+nx
for k in range (ez):
    for j in range(ey):
        for i in range(ex):
            els.append([n1,n2,n3,n4,n5,n6,n7,n8])
            n1=n1+1; n2=n2+1; n3=n3+1; n4=n4+1
            n5=n5+1; n6=n6+1; n7=n7+1; n8=n8+1
        n1=n1+1; n2=n2+1; n3=n3+1; n4=n4+1
        n5=n5+1; n6=n6+1; n7=n7+1; n8=n8+1
    n1=n1+nx; n2=n2+nx; n3=n3+nx; n4=n4+nx
    n5=n5+nx; n6=n6+nx; n7=n7+nx; n8=n8+nx

out["nodes"] = nds
out["elements"] = els
out["gauss"] = gss
out["material"] = mat
out["boundary"] = bds
out["load"] = lds

json_object = json.dumps(out, indent=4)
 
with open("sample.json", "w") as outfile:
    outfile.write(json_object)
