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
    peri_rate = 10**(x[16]) * 1.0 # periphery rate constant 
    pvav_diffuse = 40.0 # pvav diffusion rate
    center = 17
    
    al_diffuse = 0.04 # activating ligand diffusion rate
    ar_diffuse = 0.04 # activating receptor diffusion rate
    sfk_diffuse = 0.04 #sfk diffusion rate
    ir_diffuse = 0.04 # ir diffusion rate
    
    num_al = 10**(x[12])*0.5 # number of AL
    ran_comp = 10**(x[13])*1.0 # random component for microcluster movement
    k_influx = 10**(x[14])*0.001 # AR-AL influx
    
    kon_ar = 10**(x[0])* 0.002
    koff_ar = 0.023
    beta = 10**(x[1]) * 0.1  # beta
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
    kon_ar_vav_prime = 2 * kon_ar_vav
    k_vav_potn = 1.0
    k_ring = 1.0 
    k_c_diss = 0.023 
    k_p_diss = 0.023 
    k_V_diss = 0.023 
    k_CC_diss = 0.023 
    k_PP_diss = 0.023 
    k_R_diss = 0.023 
    k_DD_diss = 0.023 
    k_QQ_diss = 0.023 
            
    with open("run.sh",'w') as fRun:
        fRun.write("#!/bin/sh\n")
        fRun.write("#SBATCH --job-name=pso_nrad_clus\n")
        fRun.write("#SBATCH --time=8:00:00\n") 
        fRun.write("#SBATCH --cpus-per-task="+str(num_CPUs)+"\n")
        fRun.write("#SBATCH --ntasks=1\n")
        fRun.write("#SBATCH -o slurm.%j.out # STDOUT\n")
        fRun.write("#SBATCH -e slurm.%j.err # STDERR\n")
        fRun.write("\n")
        fRun.write("module load GCC/7.3.0-2.30 OpenMPI/3.1.1\n")
    
        for n in range(1):
            
            fInName = "in."+str(n)
            fRun.write("mpirun -np 1 ./spk_redsky < "+fInName+"\n")
        
            with open(fInName,'w') as fIn:
                fIn.write("# SPPARKS \n")
                fIn.write("\n")
                fIn.write("log      log.spparks_clus\n")
                fIn.write("seed     "+str(np.random.randint(0, 10000))+"\n")
                fIn.write("\n")
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
                fIn.write("set      i10 value 14 region global_bottom\n")
                fIn.write("set      i29 value 14 region global_bottom\n")
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
                fIn.write("add_species  r # i29 - unbounb Vav1'\n")
                fIn.write("add_species  R # i30 - bounb Vav1'\n")
                fIn.write("add_species  D # i31 - AR_p +vav'+bounb SFK complex\n")
                fIn.write("add_species  Q # i32 - phosphorylated Vav1 attached to AR_p\n")
                fIn.write("add_species  F # i33 - free/unbound phosphorylated vav'\n")
                fIn.write("add_species  Z # i34 - SHP attached to IR_p + p_Vav' attached to AR_p\n")
                fIn.write("add_species  z # i35 - SHP attached to IR_p + free p_Vav'\n")
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
                fIn.write("# AR_p + vav'\n")
                fIn.write("add_rxn      19 local p r nbr "+str(kon_ar_vav_prime)+" local R nbr\n")
                fIn.write("add_rxn      20 local R nbr "+str(koff_ar_vav)+" local p r nbr\n")
                fIn.write("\n")
                fIn.write("# ARp_vav' + SFK\n")
                fIn.write("add_rxn      21 local R k nbr "+str(kon_ar_vav_sfk)+" local D nbr\n")
                fIn.write("add_rxn      22 local D nbr "+str(koff_ar_vav_sfk)+" local R k nbr\n")
                fIn.write("\n")
                fIn.write("# Vav' phosphoyntion attached to AP_p\n")
                fIn.write("add_rxn      23 local D nbr "+str(k_ar_vav_p)+" local Q k nbr\n")
                fIn.write("\n")
                fIn.write("# dephospho pVav' attached to AP_p\n")
                fIn.write("add_rxn      24 local Q nbr "+str(k_ar_vav_dp)+" local R nbr\n")
                fIn.write("\n")
                fIn.write("# dissociation complex-phosphoVav' to AR_p + phosp-vav'\n")
                fIn.write("add_rxn      25 local Q nbr "+str(koff_ar_vav)+" local p F nbr\n")
                fIn.write("add_rxn      26 local p F nbr "+str(kon_ar_vav_prime)+" local Q nbr\n")
                fIn.write("\n")
               # fIn.write("#free pvav diffusion\n")
               # fIn.write("add_rxn      17 local f nbr b "+str(pvav_diffuse)+" local nbr f b\n")
               # fIn.write("\n")
                fIn.write("# dephospho of free pvav'\n")
                fIn.write("add_rxn      27 local F nbr "+str(k_ar_vav_dp)+" local r nbr\n")
                fIn.write("\n")
                fIn.write("# Vav dependent movement of AR (A), AR_p (p),AR-SFK (c)\n")
                fIn.write("# AR_p-vav (V;bound vav), AR_p-vav-SFK (C),AR_p-Vav_p (P)\n")
                fIn.write("add_rxn      28 local A nbr E "+str(k_vav_potn)+" local nbr A E\n")
                fIn.write("add_rxn      29 local c nbr E "+str(k_vav_potn)+" local nbr c E\n")
                fIn.write("add_rxn      30 local p nbr E "+str(k_vav_potn)+" local nbr p E\n")
                fIn.write("add_rxn      31 local V nbr E "+str(k_vav_potn)+" local nbr V E\n")
                fIn.write("add_rxn      32 local C nbr E "+str(k_vav_potn)+" local nbr C E\n")
                fIn.write("add_rxn      33 local P nbr E "+str(k_vav_potn)+" local nbr P E\n")
                fIn.write("\n")
                fIn.write("# Vav dependent movement of \n")
                fIn.write("# AR_p-vav' (R;bound vav), AR_p-vav'-SFK (D),AR_p-Vav'_p (Q)\n")
                fIn.write("add_rxn      34 local R nbr E "+str(k_vav_potn)+" local nbr R E\n")
                fIn.write("add_rxn      35 local D nbr E "+str(k_vav_potn)+" local nbr D E\n")
                fIn.write("add_rxn      36 local Q nbr E "+str(k_vav_potn)+" local nbr Q E\n")
                fIn.write("\n")
                fIn.write("# centripetal movement of bound Ar and Ar attached to Vav \n")
                fIn.write("add_rxn      37 local A nbr E "+str(k_ring)+" local nbr A E\n")
                fIn.write("add_rxn      38 local c nbr E "+str(k_ring)+" local nbr c E\n")
                fIn.write("add_rxn      39 local p nbr E "+str(k_ring)+" local nbr p E\n")
                fIn.write("add_rxn      40 local V nbr E "+str(k_ring)+" local nbr V E\n")
                fIn.write("add_rxn      41 local C nbr E "+str(k_ring)+" local nbr C E\n")
                fIn.write("add_rxn      42 local P nbr E "+str(k_ring)+" local nbr P E\n")
                fIn.write("\n")
                fIn.write("# centripetal movement of Ar attaceh to Vav' \n")
                fIn.write("add_rxn      43 local R nbr E "+str(k_ring)+" local nbr R E\n")
                fIn.write("add_rxn      44 local D nbr E "+str(k_ring)+" local nbr D E\n")
                fIn.write("add_rxn      45 local Q nbr E "+str(k_ring)+" local nbr Q E\n")
                fIn.write("\n")
                fIn.write("# centripetal movement of Vav' \n")
                fIn.write("add_rxn      46 local r nbr E "+str(k_ring)+" local nbr r E\n")
                fIn.write("add_rxn      47 local F nbr E "+str(k_ring)+" local nbr F E\n")
                fIn.write("\n")
                fIn.write("# dissociation of Ar complexes \n")
                fIn.write("add_rxn      48 local c nbr "+str(k_c_diss)+" local l a k nbr\n")
                fIn.write("add_rxn      49 local p nbr "+str(k_p_diss)+" local l a nbr\n")
                fIn.write("add_rxn      50 local V nbr "+str(k_V_diss)+" local l a v nbr\n")
                fIn.write("add_rxn      51 local C nbr "+str(k_CC_diss)+" local l a v k nbr\n")
                fIn.write("add_rxn      52 local P nbr "+str(k_PP_diss)+" local l a v nbr\n")
                fIn.write("add_rxn      53 local R nbr "+str(k_R_diss)+" local l a r nbr\n")
                fIn.write("add_rxn      54 local D nbr "+str(k_DD_diss)+" local l a r k nbr\n")
                fIn.write("add_rxn      55 local Q nbr "+str(k_QQ_diss)+" local l a r nbr\n")
                fIn.write("\n")
                fIn.write("# AR(a) influx \n")
                fIn.write("add_rxn      56 local A nbr B "+str(k_influx)+" local nbr A B\n")
                fIn.write("\n")
                fIn.write("# diffusion IR \n")
                fIn.write("add_rxn      57 local i nbr b "+str(ir_diffuse)+" local nbr i b\n")
                fIn.write("\n")
                fIn.write("\n")
                fIn.write("diag_style   propensity\n")
                fIn.write("diag_style   array i1 sum i2 sum i3 sum i4 sum i5 sum i6 sum i7 sum i8 sum i9 sum i10 sum i11 sum i12 sum i13 sum i14 sum i15 sum i16 sum i17 sum i18 sum i19 sum i20 sum i21 sum i22 sum i23 sum i24 sum i25 sum i26 sum i27 sum i28 sum i29 sum i30 sum i31 sum i32 sum i33 sum i34 sum i35 sum\n")
                fIn.write("\n")
                fIn.write("stats           1.0\n")
                fIn.write("dump            1 sites 60.0 sites."+str(n)+".* id i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35\n") 
                fIn.write("\n")
                fIn.write("run             61.0\n")
        fRun.write("sleep 30 \n")

