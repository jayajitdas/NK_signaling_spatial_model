/* -------------------------i---------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include <iostream>
#include <algorithm>

#include "math.h"
#include "string.h"
#include "app_inb_diff_custom.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{LOCAL,NBR};

/* ---------------------------------------------------------------------- */

AppInbDiffCustom::AppInbDiffCustom(SPPARKS *spk, int narg, char **arg) : 
  AppInbDiff(spk,narg,arg)
{
  // parse arguments for InbDiffCustom class only, not children
  if (strcmp(style,"inb/diff/custom") != 0) return;
  if (narg != 6) error->all(FLERR,"Illegal app_style command");
  total_cluster = atof(arg[1]);
  beta = atof(arg[2]);
  radius = atoi(arg[3]);
  ran_com = atof(arg[4]);
  periphery_rate = atof(arg[5]);
}

/* ----------------------------------------------------------------------
 * Custom rate multiplier function
------------------------------------------------------------------------- */

double AppInbDiffCustom::custom_multiplier(int i, int rstyle, int which, int jpartner)
{
  
   if (which==2) { // Diffusion: Activating ligands
     if(population[26][jpartner] != 0) {
       int ar_ligand_fill_ahead = population[0][jpartner]+ population[17][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner];
       if (ar_ligand_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
       
   } else if (which==3) { // Difusion: Activating receptors
     if(population[26][jpartner] != 0) {
       int u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[8][jpartner];
       int pot_u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[24][jpartner];
       if (u_act_fill_ahead < 2500) {
         return 1.0;
       } else return 0.0;
     } else return 0.0;
       
   } else if (which==4) { //Difussion: Src Kinase
     if(population[26][jpartner] != 0) {
      int ir_inb = population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i];
       int u_sfk_fill_ahead = population[8][jpartner]+population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner];
         if (u_sfk_fill_ahead >=2500) {
           return 0.0;
         } else {
            return 1.0;
        }
     } else return 0.0;
       
  } else if (which==5) { // Activating receptors + SFK
     if(population[1][i] != 0) {
         return 1.0;
     } else return 0.0;
   } else if (which==17) { // Diffusion:free  pVav (cytosol)
     if(population[26][jpartner] != 0) {
        return 1.0;
     } else return 0.0;
       
   } else if (which==36) { // Influx: Activating receptors
     if(population[16][i] != 0 && population[27][jpartner] != 0) {
       int u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[8][jpartner];
       if (u_act_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
       
   } else if (which==37) { // Diffusion: unbound Inhibitory Receptors
     if(population[26][jpartner] != 0) {
       int micro_inb = population[19][i] +population[20][i]+population[22][i]+population[23][i]+population[24][i]+population[25][i];
       int micro_inb_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[8][jpartner];
         if( micro_inb_fill_ahead >=2500) {
           return 0.0;
         } else return 1.0;
     } else return 0.0;
       
   } else if (which==38) { // Diffusion: Inhibitory ligands
     if(population[26][jpartner] != 0) {
       int ir_ligand_fill_ahead = population[0][jpartner]+population[17][jpartner] + population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner];
;
       if (ir_ligand_fill_ahead >=2500) {
         return 0.0;
       } else return 1.0;
     } else return 0.0;
       
  } else if (which==41) { //  Inhibitory receptors + SFK
     if(population[1][i] != 0) {
         return 1.0;
     } else return 0.0;

   } else if (which>=19 && which<=24) { // AR microcluster formation - pVAV dependent
    if ( (population[1][jpartner] != 0) && (population[1][i] != 0 )) {
      double energy_voxel = -((population[11][i]+population[14][i]+population[24][i]+population[25][i])*(population[11][i]+population[14][i]+population[24][i]+population[25][i]));
      double energy_voxel_min = 0.0;

      double energy_nbr_voxel = -((population[11][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner])*(population[11][jpartner]+population[14][jpartner]+population[24][jpartner]+population[25][jpartner]));

      double energy_nbr_voxel_min = 0.0;
      int k;
      int jj;
      for (int kk=0; kk < numneigh[i]; kk++) {
        k = neighbor[i][kk];
        double energy_voxel_temp = -((population[11][k]+population[14][k]+population[24][k]+population[25][k])*(population[11][k]+population[14][k]+population[24][k]+population[25][k]));
        energy_voxel_min = energy_voxel_temp;

        jj = neighbor[jpartner][kk];
        double energy_nbr_voxel_temp = -((population[11][jj]+population[14][jj]+population[24][jj]+population[25][jj])*(population[11][jj]+population[14][jj]+population[24][jj]+population[25][jj]));
        energy_nbr_voxel_min = energy_nbr_voxel_temp;
 
        if (energy_voxel > energy_voxel_min) {
          energy_voxel = energy_voxel_min;
        } 
        if (energy_nbr_voxel > energy_nbr_voxel_min) {
          energy_nbr_voxel = energy_nbr_voxel_min;
        }  
      }
 
      double vav_quorum;
      double vav_dep_pot;
      double total_vav_dep_pot;
      if (energy_voxel == 0 && energy_nbr_voxel == 0) {
        vav_quorum = 0.0;
      } else if (energy_voxel !=0 || energy_nbr_voxel !=0 ) {
        vav_dep_pot = (beta*(energy_voxel-energy_nbr_voxel));
        total_vav_dep_pot = exp(vav_dep_pot)/(1.0+exp(vav_dep_pot));
        vav_quorum = std::min(1.0,total_vav_dep_pot);

      }
      int lck_pot_u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[8][jpartner];
      int pot_u_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[24][jpartner];
      if (lck_pot_u_act_fill_ahead < 2500 && pot_u_act_fill_ahead < 40) {
        return vav_quorum;
      } else return 0.0;
    } else return 0.0;
       
  } else if (which>=25 && which<=30) { //Centripetal movement of AR microclusters
      if ((population[1][jpartner] != 0) && (population[1][i] !=0 )){
        double quorum_ar = 0.0;
        int coor_i;
        int coor_j;
        coor_i = population[2][i];
        coor_j = population[3][i];

// centrial movement

        if ( jpartner == i-1) {
           quorum_ar = periphery_rate * ((ran_com* ((1+(((total_pvav)/(total_cluster+total_pvav))*((coor_i/((radius-2.0)*1.0))-1)))/4.0)) +((1-ran_com)*0.25));
        } else if ( jpartner == i+1) {
           quorum_ar = periphery_rate * ((ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_i/((radius-2.0)*1.0))-1)))/4.0)) + ((1-ran_com)*0.25));
        } else if ((jpartner - i) < -1) {
           quorum_ar = periphery_rate * ((ran_com*( (1+(((total_pvav)/(total_cluster+total_pvav))*((coor_j/((radius-2.0)*1.0))-1)))/4.0)) +((1-ran_com)*0.25));
        } else if ((jpartner - i) > 1) {
           quorum_ar = periphery_rate * ((ran_com*( (1-(((total_pvav)/(total_cluster+total_pvav))*((coor_j/((radius-2.0)*1.0))-1)))/4.0)) +((1-ran_com)*0.25));
        }

        int lck_micro_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[18][jpartner]+population[19][jpartner] +population[20][jpartner]+population[22][jpartner]+population[23][jpartner]+population[24][jpartner]+population[25][jpartner]+population[8][jpartner];
        int micro_act_fill_ahead = population[5][jpartner]+population[6][jpartner]+population[7][jpartner]+population[10][jpartner]+population[11][jpartner]+population[12][jpartner]+population[13][jpartner]+population[24][jpartner];
        

// check the AR vacancy for nbr  voxel

       if ((lck_micro_act_fill_ahead < 2500) && (micro_act_fill_ahead < 40)) {
                  return quorum_ar;
        } else {
                return 0.0;
        }
    } else return 0.0;
  } else {

  return 1.0;
}
}
