#!/usr/bin/env python

import numpy as np

def setup_run(x):
    # simulation parameters
    num_CPUs = 1
    Nx = 34
    Ny = 34
    Nz = 2
    
    # biological parameters
    num_ir_cluster = 10**(x[15])* 500 # K = centripetal force component
#    central_rate = 10**(x[18])* 0.04 # central rate constant
    peri_rate = 10**(x[16]) * 1.0 # periphery rate constant
    pvav_diffuse = 40.0 #10**(x[1])* 1.0 # pvav diffuse
    center = 17
    
    al_diffuse = 0.04 # activating ligand diffuse
    ar_diffuse = 0.04 # activating receptor diffuse
    sfk_diffuse = 0.04 #sfk diffuse
    ir_diffuse = 0.04 # ir diffuse
    
    num_al = 10**(x[12])*0.5 # number of AL
 #   width = int((10**(x[13]))*1) # width at center which do no experiene centripetal force
    ran_comp = 10**(x[13])*1.0 # random component at periphery
 #   center_ran_comp = 10**(x[15])*0.1 # random component at center
    k_influx = 10**(x[14])*0.001 # AR-AL influx
    
    kon_ar = 10**(x[0])* 0.002
    koff_ar = 0.023
    beta = 10**(x[1]) *0.1  # beta
    kon_ar_sfk = 10**(x[2]) * 7.0
    koff_ar_sfk = 10**(x[3])*0.6
    k_ar_p = 10**(x[4])*0.1
    kon_ar_vav = 10**(x[5]) * 0.004
    koff_ar_vav = 10**(x[6])*1.0
    kon_ar_vav_sfk = 10**(x[7]) * 0.05
    koff_ar_vav_sfk = 10**(x[8])*1.0
    k_ar_vav_p = 10**(x[9])*1.0
    k_ar_vav_dp = 10**(x[10])*0.1
    k_ar_dp = 10**(x[11])*0.2
    #koff_par_pvav = 10**(x[12])*0.01
    #kon_par_pvav = 10**(x[13])*0.1
    #k_free_pvav_dp = 10**(x[14])*0.01
   # kon_ar_vav_prime = 2 * kon_ar_vav
    k_vav_potn = 1.0
    k_ring = 1.0
    k_c_diss = 0.023 #10**(x[15])*0.01
    k_p_diss = 0.023 # 10**(x[16])*0.01
    k_V_diss = 0.023 # 10**(x[17])* 0.01
    k_CC_diss = 0.023 # 10**(x[18])*0.01
    k_PP_diss = 0.023 # 10**(x[19])*0.01
    k_R_diss = 0.023 # 10**(x[19])*0.01
    k_DD_diss = 0.023 # 10**(x[19])*0.01
    k_QQ_diss = 0.023 # 10**(x[19])*0.01
            
    with open("run.sh",'w') as fRun:
        fRun.write("#!/bin/sh\n")
        fRun.write("#SBATCH --mail-user=Rajdeep.Grewal@nationwidechildrens.org\n")
        fRun.write("#SBATCH --mail-type=FAIL,REQUEUE\n")
        fRun.write("#SBATCH --job-name=pso_onlyAR\n")
        fRun.write("#SBATCH --time=8:00:00\n") # 7-0:00:00
        fRun.write("##SBATCH --partition=himem\n") # 7-0:00:00
        fRun.write("#SBATCH --cpus-per-task="+str(num_CPUs)+"\n")
        fRun.write("#SBATCH --ntasks=2\n")
        fRun.write("#SBATCH -o slurm.%j.out # STDOUT\n")
        fRun.write("#SBATCH -e slurm.%j.err # STDERR\n")
        fRun.write("\n")
        fRun.write("module load GCC/7.3.0-2.30 OpenMPI/3.1.1\n")
    
        for n in range(1):
            
            fInName = "in."+str(n)
            #fRun.write("mpiexec -np "+str(num_CPUs)+" ./spk_redsky < "+fInName+"\n")
            fRun.write("mpirun -np 1 ./spk_redsky < "+fInName+"\n")
        
            with open(fInName,'w') as fIn:
                fIn.write("# SPPARKS \n")
                fIn.write("\n")
                fIn.write("seed     "+str(np.random.randint(0, 10000))+"\n")
                fIn.write("\n")
                #fIn.write("app_style    inb/diff/custom "+str(num_ir_cluster)+" "+str(beta)+" "+str(center)+" "+str(width)+" "+str(ran_comp)+" "+str(center_ran_comp)+" "+str(peri_rate)+" "+str(central_rate)+"\n")
                fIn.write("app_style    inb/diff/custom "+str(num_ir_cluster)+" "+str(beta)+" "+str(center)+" "+str(ran_comp)+" "+str(peri_rate)+"\n")
                fIn.write("\n")
                fIn.write("dimension    3\n")
                fIn.write("lattice      sc/6n 0.5\n")
                fIn.write("region       box block 0 "+str(Nx)+" 0 "+str(Ny)+" 0 "+str(Nz)+"\n")
                fIn.write("create_box   box\n")
                fIn.write("create_sites box\n")
                fIn.write("\n")
                fIn.write("region       sim_box block 2 31 2 31 0 0\n")
                fIn.write("region       no_boundary_sim_box block 3 30 3 30 0 0\n")
                fIn.write("region       global_bottom block 0 34 0 34 0 0\n")
                fIn.write("region       cytosol block 2 31 2 31 0 0\n")
                fIn.write("region       boundary_wall intersect 2 global_bottom no_boundary_sim_box side out\n")
                fIn.write("region       boundary_sim_box intersect 2 boundary_wall sim_box side in\n")
                fIn.write("region       wall intersect 2 global_bottom sim_box side out\n")
                fIn.write("region       cyt_box block 0 34 0 34 1 1\n")
                fIn.write("\n")
                fIn.write("read_sites       sites.30.30\n")
                fIn.write("set      i1 value 3 region global_bottom fraction "+str(num_al)+"\n")
                fIn.write("set      i2 value 1 region sim_box\n")
                fIn.write("set      i6 value 2 region global_bottom\n")
                fIn.write("set      i9 value 175 region global_bottom\n")
                fIn.write("set      i10 value 28 region global_bottom\n")
                #fIn.write("set      i29 value 14 region global_bottom\n")
                fIn.write("set      i17 value 1 region wall\n")
                fIn.write("set      i17 value 0 region cyt_box\n")
                fIn.write("set      i19 value 26 region global_bottom\n")
                fIn.write("set      i22 value 522 region global_bottom\n")
                fIn.write("set      i27 value 1 region global_bottom\n")
                fIn.write("set      i28 value 1 region boundary_sim_box\n")
                fIn.write("\n")
                fIn.write("sector   yes\n")
                fIn.write("solve_style  tree \n")
                fIn.write("\n")
                fIn.write("# species (note that in the C++ code indexing starts at 0 instead of 1)\n")
                fIn.write("add_species  l # i1 - unbound activating ligand\n")
                fIn.write("add_species  E # i2 - empty space (i.e. not a wall or free surface)\n")
                fIn.write("add_species  y # i3 - X cord \n")
                fIn.write("add_species  Y # i4 - Y cord \n")
                fIn.write("add_species  N # i5 - Voxel No.a\n")
                fIn.write("add_species  a # i6 - unbound activating receptors\n")
                fIn.write("add_species  A # i7 - bound activating receptors\n")
                fIn.write("add_species  p # i8 - phosphorylated bound activating receptors\n")
                fIn.write("add_species  k # i9 - SFK kinases\n")
                fIn.write("add_species  v # i10 -  unbounb Vav1s\n")
                fIn.write("add_species  V # 11 -  bounb Vav1s\n")
                fIn.write("add_species  P # i12 - phosphorylated Vav1 attached to AR_p\n")
                fIn.write("add_species  c # i13 - Lck+A complex \n")
                fIn.write("add_species  C # i14 - AR_p +vav+bounb SFK complex \n")
                fIn.write("add_species  f # i15 - free/unbound phosphorylated vav\n")
                fIn.write("add_species  W # i16 - wall\n")
                fIn.write("add_species  w # i17 - influx wall\n")
                fIn.write("add_species  L # i18 - unbound inhibiting ligands\n")
                fIn.write("add_species  i # i19 - unbound inhibiting receprors\n")
                fIn.write("add_species  I # i20 - bound inhibiting receprors\n")
                fIn.write("add_species  h # i21 - phosphorylated bound inhibiting receprors\n")
                fIn.write("add_species  s # i22 - unbound SHP\n")
                fIn.write("add_species  S # i23 - bound SHP\n")
                fIn.write("add_species  x # i24 - Lck+I complex\n")
                fIn.write("add_species  G # i25 - SHP attached to IR_p + p_Vav attached to AR_p\n")
                fIn.write("add_species  g # i26 - SHP attached to IR_p + free p_Vav\n")
                fIn.write("add_species  b # i27 - global sim box\n")
                fIn.write("add_species  B # i28 - boundary of sim_box (next to wall)\n")
                fIn.write("\n")
                fIn.write("# #AR + L_A\n")
                fIn.write("add_rxn      0 local l a nbr "+str(kon_ar)+" local A nbr\n")
                fIn.write("add_rxn      1 local A nbr "+str(koff_ar)+" local l a nbr\n")
                fIn.write("\n")
                fIn.write("# diffusion AR , SFK(s) vav(v)\n")
                fIn.write("add_rxn      2 local l nbr b "+str(al_diffuse)+" local nbr l b\n")
                fIn.write("add_rxn      3 local a nbr b "+str(ar_diffuse)+" local nbr a b\n")
                fIn.write("add_rxn      4 local k nbr b "+str(sfk_diffuse)+" local nbr k b\n")
                fIn.write("\n")
                fIn.write("# AR+kinase\n")
                fIn.write("add_rxn      5 local A k nbr "+str(kon_ar_sfk)+" local c nbr\n")
                fIn.write("add_rxn      6 local c nbr "+str(koff_ar_sfk)+" local A k nbr\n")
                fIn.write("\n")
                fIn.write("# phophorylation of  AR\n")
                fIn.write("add_rxn      7 local c nbr "+str(k_ar_p)+" local p k nbr\n")
                fIn.write("\n")
                fIn.write("# AR_p + vav\n")
                fIn.write("add_rxn      8 local p v nbr "+str(kon_ar_vav)+" local V nbr\n")
                fIn.write("add_rxn      9 local V nbr "+str(koff_ar_vav)+" local p v nbr\n")
                fIn.write("\n")
                fIn.write("# ARp_vav + SFK\n")
                fIn.write("add_rxn      10 local V k nbr "+str(kon_ar_vav_sfk)+" local C nbr\n")
                fIn.write("add_rxn      11 local C nbr "+str(koff_ar_vav_sfk)+" local V k nbr\n")
                fIn.write("\n")
                fIn.write("# Vav phosphoyntion attached to AP_p\n")
                fIn.write("add_rxn      12 local C nbr "+str(k_ar_vav_p)+" local P k nbr\n")
                fIn.write("\n")
                fIn.write("# dephospho\n")
                fIn.write("add_rxn      13 local P nbr "+str(k_ar_vav_dp)+" local V nbr\n")
                fIn.write("add_rxn      14 local p nbr "+str(k_ar_dp)+" local A nbr\n")
                fIn.write("\n")
                fIn.write("# dissociation complex-phosphoVav to AR_p + phosp-vav\n")
                fIn.write("add_rxn      15 local P nbr "+str(koff_ar_vav)+" local p f nbr\n")
                fIn.write("add_rxn      16 local p f nbr "+str(kon_ar_vav)+" local P nbr\n")
                fIn.write("\n")
                fIn.write("#free pvav diffusion\n")
                fIn.write("add_rxn      17 local f nbr b "+str(pvav_diffuse)+" local nbr f b\n")
                fIn.write("\n")
                fIn.write("# dephospho of free pvav\n")
                fIn.write("add_rxn      18 local f nbr "+str(k_ar_vav_dp)+" local v nbr\n")
                fIn.write("\n")
                fIn.write("# Vav dependent movement of AR (A), AR_p (p),AR-SFK (c)\n")
                fIn.write("# AR_p-vav (V;bound vav), AR_p-vav-SFK (C),AR_p-Vav_p (P)\n")
                fIn.write("add_rxn      19 local A nbr E "+str(k_vav_potn)+" local nbr A E\n")
                fIn.write("add_rxn      20 local c nbr E "+str(k_vav_potn)+" local nbr c E\n")
                fIn.write("add_rxn      21 local p nbr E "+str(k_vav_potn)+" local nbr p E\n")
                fIn.write("add_rxn      22 local V nbr E "+str(k_vav_potn)+" local nbr V E\n")
                fIn.write("add_rxn      23 local C nbr E "+str(k_vav_potn)+" local nbr C E\n")
                fIn.write("add_rxn      24 local P nbr E "+str(k_vav_potn)+" local nbr P E\n")
                fIn.write("\n")
                fIn.write("# centripetal movement of bound Ar and Ar attached to Vav \n")
                fIn.write("add_rxn      25 local A nbr E "+str(k_ring)+" local nbr A E\n")
                fIn.write("add_rxn      26 local c nbr E "+str(k_ring)+" local nbr c E\n")
                fIn.write("add_rxn      27 local p nbr E "+str(k_ring)+" local nbr p E\n")
                fIn.write("add_rxn      28 local V nbr E "+str(k_ring)+" local nbr V E\n")
                fIn.write("add_rxn      29 local C nbr E "+str(k_ring)+" local nbr C E\n")
                fIn.write("add_rxn      30 local P nbr E "+str(k_ring)+" local nbr P E\n")
                fIn.write("\n")
                fIn.write("# dissociation of Ar complexes \n")
                fIn.write("add_rxn      31 local c nbr "+str(k_c_diss)+" local l a k nbr\n")
                fIn.write("add_rxn      32 local p nbr "+str(k_p_diss)+" local l a nbr\n")
                fIn.write("add_rxn      33 local V nbr "+str(k_V_diss)+" local l a v nbr\n")
                fIn.write("add_rxn      34 local C nbr "+str(k_CC_diss)+" local l a v k nbr\n")
                fIn.write("add_rxn      35 local P nbr "+str(k_PP_diss)+" local l a v nbr\n")
                fIn.write("\n")
                fIn.write("# AR(a) influx \n")
                fIn.write("add_rxn      36 local A nbr B "+str(k_influx)+" local nbr A B\n")
                fIn.write("\n")
                #fIn.write("# IR(i) + I_L (L) \n")
                #fIn.write("add_rxn      38 local L i nbr 1.0 local I nbr\n")
                #fIn.write("add_rxn      39 local I nbr 0.01 local L i nbr\n")
                #fIn.write("\n")
                fIn.write("# IR(i) + I_L (L) \n")
                #fIn.write("add_rxn      40 local L nbr E 0.01 local nbr L E\n")
                fIn.write("add_rxn      37 local i nbr b "+str(ir_diffuse)+" local nbr i b\n")
                #fIn.write("#add_rxn      39 local s nbr E "+str(shp_diffuse)+" local nbr s E\n")
                fIn.write("\n")
                fIn.write("\n")
                fIn.write("diag_style   propensity\n")
                fIn.write("diag_style   array i1 sum i2 sum i3 sum i4 sum i5 sum i6 sum i7 sum i8 sum i9 sum i10 sum i11 sum i12 sum i13 sum i14 sum i15 sum i16 sum i17 sum i18 sum i19 sum i20 sum i21 sum i22 sum i23 sum i24 sum i25 sum i26 sum i27 sum i28 sum\n")
                fIn.write("\n")
                fIn.write("stats           1.0\n")
                fIn.write("dump            1 sites 60.0 sites."+str(n)+".* id i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 i21 i22 i23 i24 i25 i26 i27 i28\n") 
                fIn.write("\n")
                fIn.write("run             62.0\n")
        fRun.write("sleep 30 \n")
        fRun.write("if [ ! -f sites.0.1 ]; then echo \"File not created \" 1>&2 ; fi \n")
        #fRun.write("if [ ! -f sites.0.3 ]; then echo \"File not created \" 1>&2 ; fi \n")
        #fRun.write("if [ ! -f sites.0.5 ]; then echo \"File not created \" 1>&2 ; fi \n")
        #fRun.write("if [ ! -f sites.0.7 ]; then echo \"File not created \" 1>&2 ; fi \n")
