#!/usr/bin/env python3

import sys
import numpy as np

########## read in bacteria data from stdin #########

# header
for n in range(5):
        input()

line = input()
split_line = line.split()
Nx = int(split_line[1])*2

line = input()
split_line = line.split()
Ny = int(split_line[1])*2

line = input()
split_line = line.split()
Nz = int(float(split_line[1])*2)

Vol = (Nx-4)*(Ny-4)*Nz

for n in range(3):
        input()

# data
nkg2d = np.zeros((Nx-4,Ny-4,Nz))

for n in range(Nx*Ny*Nz):
        line = input()
        split_line = line.split()
        site_id = int(split_line[0]) - 1
        i = site_id%Nx
        j = site_id//Nx%Ny
        k = site_id//(Nx*Ny)
        if (i != 0 and i != 1 and j != 1 and j != 0 and i != Nx-2 and i != Nx-1 and j != Ny-2 and j != Ny-1):
                nkg2d[int(split_line[3])-1,int(split_line[4])-1,k] = int(split_line[6])+ int(split_line[7])  + int(split_line[8]) + int(split_line[11]) + int(split_line[12]) + int(split_line[13]) + int(split_line[14])+ int(split_line[25])+ int(split_line[30])+int(split_line[31])+int(split_line[32])+int(split_line[34])
 

########## print in vtk format to stdout #########
sys.stdout.write("# vtk DataFile Version 2.0\n");
sys.stdout.write("InbSimInVitro\n");
sys.stdout.write("ASCII\n");
sys.stdout.write("DATASET STRUCTURED_POINTS\n");
sys.stdout.write("DIMENSIONS %d %d %d\n" % (Nx-4,Ny-4,Nz));
sys.stdout.write("ASPECT_RATIO 1 1 1\n");
sys.stdout.write("ORIGIN 0 0 0\n");
sys.stdout.write("POINT_DATA %d\n" % Vol);
sys.stdout.write("SCALARS u double 1\n");
sys.stdout.write("LOOKUP_TABLE default\n");
# VTK file format is in column major order
for k in range(Nz):
        for j in range(Ny-4):
                for i in range(Nx-4):
                        sys.stdout.write("%f\n" % bact_ctr[i,j,k])



###### command to run the file

# for i in `seq 0 180`;do if [ -f sites.0.$i ];then echo "sites.0.$i => ar_only.$i.vtk"; if [ ! -f ar_only.$i.vtk ]; then ../sites2vtk_all_ar.py <sites.0.$i > ar_only.$i.vtk; fi; fi;done


### The above command will convert each sites.*.* output fils from the model runs to corresponding ar_only.*.vtk.


