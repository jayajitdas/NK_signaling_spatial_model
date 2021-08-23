#!/usr/bin/env python3
import os, sys, shutil
#if sys.version_info[0]<3:
#   sys.stderr.write("You need python 3 or later to run this script\n")
#   exit(1)

import numpy as np

def adjust_voxel_size(f,dirname,thres,ti):
    ########## read in NKG2D data from file #########
    
    # header
    for n in range(5):
        f.readline()
    
    line = f.readline()
    split_line = line.split()
    Nx = int(split_line[1])*2
    
    line = f.readline()
    split_line = line.split()
    Ny = int(split_line[1])*2
    
    line = f.readline()
    split_line = line.split()
    Nz = int(split_line[1])*2
    
    
    for n in range(3):
        f.readline()
    
    # data
    nk_ar = []
    id_site = np.zeros(Nx*Ny)
    val_site = []
    ll = 0
    for n in range(Nx*Ny*Nz):
        line = f.readline()
        split_line = line.split()
        site_id = int(split_line[0]) - 1
        i = site_id%Nx
        j = (site_id//Nx)%Ny
        k = site_id//(Nx*Ny)
        if (k == 0 and i != 0 and i != 1 and j != 1 and j != 0 and i != Nx-2 and i != Nx-1 and j != Ny-2 and j != Ny-1): # simulation box of size 15 micron x 15 micron
            
            nk_ar.append( float(split_line[6]) + float(split_line[7]) + int(split_line[8]) + int(split_line[11]) + int(split_line[12]) + int(split_line[13]) + int(split_line[14])+ int(split_line[30])+int(split_line[31])+int(split_line[32]))
            ll = ll+1
                
    f1 = open("sites."+str(thres)+"_"+str(dirname)+"adjusted15.15"+"_"+str(ti),"w")
    f1.write("Site file written by dump sites 1 command\n\n")
    f1.write("3 dimension\n")
    f1.write(str(int((Ny-4)*(Nx-4)))+" sites\n")
    f1.write("id i1 values\n")
    f1.write("0 "+str(Nx-4)+" xlo xhi\n")
    f1.write("0 "+str(Ny-4)+" ylo yhi\n")
    f1.write("0 "+str(Nz)+" zlo zhi\n\n")
    f1.write("Values\n\n")
    term = 0.0
      
    for ii in range(len(val_site)):
        if(nk_ar[ii] >= thres):
            f1.write("%d %f\n"%(ii+1, nk_ar[ii]))
        else:
            f1.write("%d 0.0\n"%(ii+1))
   
    f1.close()
    if os.path.isfile(dirname+"/sites."+str(thres)+"_"+str(dirname)+"adjusted15.15"+"_"+str(ti)):
        os.remove("sites."+str(thres)+"_"+str(dirname)+"adjusted15.15"+"_"+str(ti))
    else:
        shutil.move("sites."+str(thres)+"_"+str(dirname)+"adjusted15.15"+"_"+str(ti),dirname)

    return (None)
