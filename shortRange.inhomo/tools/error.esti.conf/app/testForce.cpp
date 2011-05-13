#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>

#include "ForceKernel.h"

int main(int argc, char * argv[])
{
  double dxv = 0.1;
  if (argc != 1){
    dxv = atof(argv[1]);
  }
  
  double ref[3];
  double dx[3];
  double target0[3];
  double target1[3];
  for (int i = 0; i < 3; ++i) {
    ref[i] = 5.;
    dx[i] = dxv;
  }
  dx[1] = dx[2] = 0.;
  // dx[0] = -dxv;
  for (int i = 0; i < 3; ++i) {
    target0[i] = ref[i] + dx[i];
    target1[i] = ref[i] - dx[i];
  }
  double ref_force[3];
  double target_force0[3];
  double target_force_taylor0[3];
  double target_force1[3];
  double target_force_taylor1[3];

  Disperson6_Taylor force (1., 1., 1.);

  force.f(ref[0], ref[1], ref[2],
	  ref_force[0],
	  ref_force[1],
	  ref_force[2]);
  force.f (target0[0], target0[1], target0[2],
	   target_force0[0],
	   target_force0[1],
	   target_force0[2]);
  force.f (target1[0], target1[1], target1[2],
	   target_force1[0],
	   target_force1[1],
	   target_force1[2]);  
  force.f_Taylor (ref[0], ref[1], ref[2],
		  dx[0], dx[1], dx[2],
		  target_force_taylor0[0],
		  target_force_taylor0[1],
		  target_force_taylor0[2]);
  force.f_Taylor (ref[0], ref[1], ref[2],
		  -dx[0], -dx[1], -dx[2],
		  target_force_taylor1[0],
		  target_force_taylor1[1],
		  target_force_taylor1[2]);
  // test group 1:
  printf ("////////////////////////////////////////////////////////////\n");
  {
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = force.fdotf_oneSide_taylor (ref[0], ref[1], ref[2],
					      dx[0], dx[1], dx[2]);
    for (int i = 0; i < 3; ++i){
      sum0 += target_force0[i] * ref_force[i];
      sum1 += target_force_taylor0[i] * ref_force[i];
    }
    printf ("ref:     %e %e %e\n", ref[0], ref[1], ref[2]);
    printf ("target0: %e %e %e\n", target0[0], target0[1], target0[2]);
    printf ("reff     %e %e %e\n", ref_force[0], ref_force[1], ref_force[2]);
    printf ("f0 exact: (%e %e %e), taylor (%e %e %e)\n",
	    target_force0[0], target_force0[1], target_force0[2], 
	    target_force_taylor0[0], target_force_taylor0[1], target_force_taylor0[2]);
    printf ("exact: (%e), taylor1 (%e), error: %e\n",
	    sum0, sum1, fabs(sum0 - sum1));
    printf ("exact: (%e), taylor2 (%e), error: %e\n",
	    sum0, sum2, fabs(sum0 - sum2));
    printf ("\n\n");
  }
  
  // test group 2:
  {
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = force.fdotf_oneSide_taylor (ref[0], ref[1], ref[2],
					      -dx[0], -dx[1], -dx[2]);
    for (int i = 0; i < 3; ++i){
      sum0 += target_force1[i] * ref_force[i];
      sum1 += target_force_taylor1[i] * ref_force[i];
    }
    printf ("ref:     %e %e %e\n", ref[0], ref[1], ref[2]);
    printf ("target1: %e %e %e\n", target1[0], target1[1], target1[2]);
    printf ("reff     %e %e %e\n", ref_force[0], ref_force[1], ref_force[2]);
    printf ("f1 exact: (%e %e %e), taylor (%e %e %e)\n",
	    target_force1[0], target_force1[1], target_force1[2], 
	    target_force_taylor1[0], target_force_taylor1[1], target_force_taylor1[2]);
    printf ("exact: (%e), taylor1 (%e), error: %e\n",
	    sum0, sum1, fabs(sum0 - sum1));
    printf ("exact: (%e), taylor2 (%e), error: %e\n",
	    sum0, sum2, fabs(sum0 - sum2));
    printf ("\n\n");
  }

  
  {
    double sum0 = 0.;
    double sum1 = 0.;
    double sum2 = force.fdotf_taylor (ref[0], ref[1], ref[2],
				      dx[0], dx[1], dx[2]);
    for (int i = 0; i < 3; ++i){
      // double tmp = target_force_taylor[0] -  target_force[0];
      sum0 += target_force0[i] * target_force1[i];
      sum1 += target_force_taylor0[i] * target_force_taylor1[i];
    }
    sum0 = (sum0);
    sum1 = (sum1);

    printf ("r0: %e %e %e\n", target0[0], target0[1], target0[2]);
    printf ("r1: %e %e %e\n", target1[0], target1[1], target1[2]);
    printf ("f0 exact: (%e %e %e), taylor (%e %e %e)\n",
	    target_force0[0], target_force0[1], target_force0[2], 
	    target_force_taylor0[0], target_force_taylor0[1], target_force_taylor0[2]);
    printf ("f1 exact: (%e %e %e), taylor (%e %e %e)\n",
	    target_force1[0], target_force1[1], target_force1[2], 
	    target_force_taylor1[0], target_force_taylor1[1], target_force_taylor1[2]);
  
    printf ("exact: (%e), taylor1 (%e), error: %e\n",
	    sum0, sum1, fabs(sum0 - sum1));
    printf ("exact: (%e), taylor2 (%e), error: %e\n",
	    sum0, sum2, fabs(sum0 - sum2));
  }
  printf ("////////////////////////////////////////////////////////////\n");
  
  return 0;
}



// int main(int argc, char * argv[])
// {
//   double dxv = 0.1;
//   if (argc != 1){
//     dxv = atof(argv[1]);
//   }
  
//   double ref[3];
//   double dx[3];
//   double target[3];
//   for (int i = 0; i < 3; ++i) {
//     ref[i] = 5.;
//     dx[i] = dxv;
//   }
//   for (int i = 0; i < 3; ++i) {
//     target[i] = ref[i] + dx[i];
//   }
//   double target_force[3];
//   double target_force_taylor[3];

//   Disperson6_Taylor force (1., 1., 3.);
  
//   force.f (target[0], target[1], target[2],
// 	   target_force[0],
// 	   target_force[1],
// 	   target_force[2]);
//   force.f_Taylor (ref[0], ref[1], ref[2],
// 		  dx[0], dx[1], dx[2],
// 		  target_force_taylor[0],
// 		  target_force_taylor[1],
// 		  target_force_taylor[2]);

//   double sum = 0.;
//   double sum1 = 0.;
//   for (int i = 0; i < 3; ++i){
//     double tmp = target_force_taylor[i] -  target_force[i];
//     sum += tmp * tmp;
//     sum1 += target_force[i] * target_force[i];
//   }
//   sum = sqrt(sum);
//   sum1 = sqrt(sum1);

//   printf ("exact: (%e %e %e), taylor (%e %e %e), error: %e, forceMagnitude: %e\n",
// 	  target_force[0], target_force[1], target_force[2],
// 	  target_force_taylor[0], target_force_taylor[1], target_force_taylor[2],
// 	  sum, sum1);

//   return 0;
// }

