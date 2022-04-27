#!/usr/bin/env python3
import sys
#if sys.version_info[0]<3:
#   sys.stderr.write("You need python 3 or later to run this script\n")
#   exit(1)

import numpy as np

def calc_pair_corr(f):
    ########## read in NKG2D data from file #########
    
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
    nk_ar = np.zeros((Nx,Ny))
    for n in range(Nx*Ny):
        line = f.readline()
        split_line = line.split()
        site_id = int(split_line[0]) - 1
        i = site_id%Nx
        j = site_id//Nx%Ny
        k = site_id//(Nx*Ny)
        if (k == 0):
            nk_ar[i,j] = float(split_line[1])    
    
    ########## Calculate C(r,t) ##########
    # set up some helper functions
    def periodic_dist(pt1,pt2):
        dxSq = min( (pt1[0]-pt2[0])**2, (pt1[0]-pt2[0]-Nx)**2, (pt1[0]-pt2[0]+Nx)**2 )
        dySq = min( (pt1[1]-pt2[1])**2, (pt1[1]-pt2[1]-Ny)**2, (pt1[1]-pt2[1]+Ny)**2 )
        return np.sqrt(dxSq+dySq)
    
    # calculate C(r,t)
    corr_nk_ar = {}
    count_bin = {}
    rMax = min(Nx,Ny)/2
    
    # subtract off means and divide out standard deviations
    nk_ar = nk_ar - nk_ar.mean()
    nk_ar = nk_ar / nk_ar.std()
    
    for rx in range(Nx):
        for ry in range(Ny):
            dist = periodic_dist([0,0],[rx,ry])
            rdist = np.round(dist,decimals=6)
            if rdist <= rMax:
                nk_ar_overlap = np.mean( nk_ar * np.roll(nk_ar,(rx,ry),axis=(0,1)) )
                   
                if rdist not in count_bin:
                    count_bin[rdist] = 0
                    corr_nk_ar[rdist] = 0.
 
                count_bin[ rdist ] += 1
                corr_nk_ar[ rdist ] += nk_ar_overlap
    
    for key in sorted(count_bin):
        corr_nk_ar[key] /= count_bin[key]

    r = [] 
    C_nk_ar = []
    i = 0
    for key in sorted(count_bin):
        r.append(key)
        C_nk_ar.append(corr_nk_ar[key])
    
    return (C_nk_ar)

