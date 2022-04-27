#!/usr/bin/env python3

import numpy as np

def calc_S4_std_density(f):
    for n in range(5):
        f.readline()
    
    line = f.readline()
    split_line = line.split()
    Nx = int(split_line[1])
    
    line = f.readline()
    split_line = line.split()
    Ny = int(split_line[1])
    
    line = f.readline()
    split_line = line.split()
    Nz = int(split_line[1])
    
    data_dap10 = np.zeros((Nx,Ny,Nz))
    
    for n in range(3):
        f.readline()
    
    for n in range(Nx*Ny):
        line = f.readline()
        split_line = line.split()
        site_id = int(split_line[0]) - 1
        i = site_id%Nx
        j = site_id//Nx%Ny
        k = site_id//(Nx*Ny)
        data_dap10[i,j,k] = float(split_line[1])


     
    max_density = 0.0 
    for i in range(Nx):
        for j in range(Ny):
            density_dap10 = 0.0
            density_dap10 = data_dap10[i,j,0]
            if (density_dap10 > max_density):
                max_density = density_dap10
    height = np.zeros(Nx*Ny)
    ll = 0
    for i in range(Nx):
        for j in range(Ny):
            height[ll] = (( data_dap10[i,j,0] * 1.0)/(1.0 * max_density))
            ll = ll+1
    
    height_mean = height.mean()
    height_stdev = height.std()
    
    return (height_mean,height_stdev)


