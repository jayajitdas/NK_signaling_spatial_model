#!/usr/bin/env python3
import os, sys, shutil
#if sys.version_info[0]<3:
#   sys.stderr.write("You need python 3 or later to run this script\n")
#   exit(1)

import numpy as np

def S4_adjust_voxel_size(f,dirname,thres,ti):
    ########## read in bacteria data from file #########
    
    # header
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
    
    for n in range(3):
        f.readline()
    
    # data
    
    id_site = np.zeros(Nx*Ny)
    dap_site = []
    ll = 0
    for n in range(Nx*Ny):
        line = f.readline()
        split_line = line.split()
        site_id = int(split_line[0]) - 1
        i = site_id%Nx
        j = (site_id//Nx)%Ny
        k = site_id//(Nx*Ny)
        
        if (k == 0):
            dap_site.append(float(split_line[1]))
            ll = ll+1
                
    f1 = open("sites."+str(thres)+"_"+str(dirname)+"_thres_adjusted_S4_"+str(ti)+"_26.26","w")
    f1.write("Site file written by dump sites 1 command\n\n")
    f1.write("3 dimension\n")
    f1.write(str(int(((Ny)/2)*((Nx)/2)))+" sites\n")
    f1.write("id i1 values\n")
    f1.write("0 "+str(int((Nx)/2))+" xlo xhi\n")
    f1.write("0 "+str(int((Ny)/2))+" ylo yhi\n")
    f1.write("0 "+str(Nz)+" zlo zhi\n\n")
    f1.write("Values\n\n")
    term = 0.0
    
    kk = 0
    dap10 = []
    Nx1 = int((Nx)/2)
    Ny1 = int((Ny)/2)
    dap_array  = np.zeros((Nx1,Ny))
    
    # dimension reduction across x-axis 
    for ii in range(len(dap_site)):
        i = ii%(Nx)
        new_i = int(i/2)
        j = (ii//(Nx))%(Ny)
        k = ii//((Nx)*(Ny))
        if (ii <(len(dap_site))-1 ):
            if ((ii+2)%2 == 0):
                term = (dap_site[ii] + dap_site[ii+1])/2.0
                dap10.append(term)
                dap_array[new_i,j]=term
                kk = kk+1
    
    # dimension reduction across y-axis
    kkk = 0
    fdap_array  = np.zeros((Nx1,Ny1))
    dap_1 = []
    for i in range(Nx1):
        coun = 0
        for j in range(Ny):
            jj = (i//Nx1)%Ny
            new_j = int(jj/2)
            if (j < Ny-1 ):
                if ((j+2)%2 == 0):
                    term = (dap_array[i,j] + dap_array[i,j+1])/2.0
                    dap_1.append(term)
                    term2 = dap_1[kkk]/1.0
                    
                    fdap_array[i,coun] = dap_1[kkk]
                    kkk = kkk+1
                    coun = coun+1
    k2 = 0
    for j in range(Ny1):
        for i in range(Nx1):
            f1.write("%d %f\n"%(k2+1, fdap_array[i,j]))
            k2 = k2+1
    f1.close()
    
    return (None)

