#!/usr/bin/env python3

import numpy as np

def calc_std_density(f):
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
    
    nk_ar_data = np.zeros((Nx,Ny,Nz))
    
    for n in range(3):
        f.readline()
    
    for n in range(Nx*Ny):
        line = f.readline()
        split_line = line.split()
        site_id = int(split_line[0]) - 1
        i = site_id%Nx
        j = site_id//Nx%Ny
        k = site_id//(Nx*Ny)
        
        nk_ar_data[i,j,k] = float(split_line[1]) 


     
    max_density = 0.0
    for i in range(Nx):
        for j in range(Ny):
            nk_ar_density = 0.0
            nk_ar_density = nk_ar_data[i,j,0]
            if (nk_ar_density > max_density):
                max_density = nk_ar_density
    height = np.zeros(Nx*Ny)
    ll = 0
    for i in range(Nx):
        for j in range(Ny):
            height[ll] = (( nk_ar_data[i,j,0] * 1.0)/(1.0 * max_density))
            ll = ll+1
    
    height_mean = height.mean()
    height_stdev = height.std()
    
    return (height_mean,height_stdev)


