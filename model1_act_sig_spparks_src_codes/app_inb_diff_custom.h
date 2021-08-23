/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(inb/diff/custom,AppInbDiffCustom)

#else

#ifndef SPK_APP_INB_DIFF_CUSTOM_H
#define SPK_APP_INB_DIFF_CUSTOM_H

#include "app_inb_diff.h"

namespace SPPARKS_NS {

class AppInbDiffCustom : public AppInbDiff {
 public:
  AppInbDiffCustom(class SPPARKS *, int, char **);
  ~AppInbDiffCustom() {}

 protected:
  double total_cluster;
  //int food_left;
  double beta;
  int radius;
  //int center_vicinity;
  //float frac_vel;
  float ran_com;
//  float ran_com1;
  float periphery_rate;
//  float center_rate;
  //int bir_conc;
  //int *ir_cluster_center;
  double custom_multiplier(int, int, int, int);
  //void cluster_movement(int, int, int, int);
  

 // void stats(int *ir_cluster_center);
/*  void stats_header(char *);*/
};

}

#endif
#endif
