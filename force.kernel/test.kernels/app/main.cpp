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

using namespace std;

#include "nb_kernels.h"
#include "nb_interactions.h"
#include "RandomGenerator.h"
#include "Stopwatch.h"

const int MDDIM = 3;
typedef float ValueType;

int main(int argc, char * argv[])
{
  int n0 = 10000;
  int n1 = 4*10000;
  ValueType cutoff[1] = {5.0};
  ValueType param [4] = {0.};
  param[0] = 24;
  param[1] = 48;

  ValueType * coord = (ValueType *) malloc (sizeof(ValueType) * MDDIM * n0);
  int * vdwtype = (int *) malloc (sizeof(ValueType) * MDDIM * n0);
  for (int ii = 0; ii < n0 * MDDIM; ++ii){
    coord[ii] = RandomGenerator_MT19937::genrand_real3();
    vdwtype[ii] = 0;
  }
  ValueType * force = (ValueType *) malloc (sizeof(ValueType) * MDDIM * n0);
  ValueType * force_shift = (ValueType *) malloc (sizeof(ValueType) * MDDIM * n0);
  for (int ii = 0; ii < n0 * MDDIM; ++ii){
    force[ii] = force_shift[ii] = 0;
  }
  int * neighbor_index = (int *) malloc (sizeof(int) * n1);
  for (int ii = 0; ii < n1; ++ii){
    neighbor_index[ii] = int(RandomGenerator_MT19937::genrand_real2() * n0);
  }
  int * neighbor_index_data = (int *) malloc (sizeof(int) * n0 * 2);
  for (int ii = 0; ii < n0 * 2; ii+=2){
    neighbor_index_data[ii+0] = 0;
    neighbor_index_data[ii+1] = n1;
  }

  Stopwatch mywatch;

  mywatch.start();
  nb_pair_force <3,
		 float,
		 nb_interaction_geometric_none_tag,
		 nb_interaction_accelteration_128s_tag,
		 nb_interaction_electrostatic_null_tag,
		 nb_interaction_vanderWaals_cutoff_tag,
		 nb_interaction_compute_f_tag>
      (NULL, param, cutoff, n0, neighbor_index_data, neighbor_index, coord, NULL, vdwtype, 1, force, force_shift, NULL);
  mywatch.stop();
  
  for (int ii = 0; ii < n0; ++ii){
    printf ("%f %f %f\n", force[0+3*ii], force[1+3*ii], force[2+3*ii]);
  }
  
  printf ("syste: %f  user: %f  real: %f\n", mywatch.system(), mywatch.user(), mywatch.real());

  free (coord);
  free (neighbor_index);
  free (neighbor_index_data);
  free (force_shift);
  free (force);

  return 0;
}


// int main(int argc, char * argv[])
// {
//   double eleParam[5] = {0.};
//   double vdwParam[LennardJones6_12::nparam];
//   vdwParam[LennardJones6_12::epsilon] = 1.0;
//   vdwParam[LennardJones6_12::sigma] = 1.0;
  
//   double cutoff[1] = {5.0};
//   double dof[MDDIM*2];
//   dof[0] = 1.5;
//   dof[1] = 0.0;
//   dof[2] = 0.0;
//   dof[3] = 0.0;
//   dof[4] = 0.0;
//   dof[5] = 0.0;
//   int index0 = 0;
//   int index1 = 1;
//   double charge[2];
//   charge[0] = -1;
//   charge[1] = -1;  
//   double force[MDDIM*2];
//   force[0] = force[1] = force[2] = force[3] = force[4] = force[5] = 0.;

//   int n0 = 1000;
//   int n1 = 100000;

//   for (int ii = 0; ii < n0; ++ii){
//     for (int jj = 0; jj < n1; ++jj){
//       dof[0] = RandomGenerator_MT19937::genrand_real3();
// // dof1[1] = RandomGenerator_MT19937::genrand_real3();
// // dof1[2] = RandomGenerator_MT19937::genrand_real3();
// // dof0[0] = RandomGenerator_MT19937::genrand_real3();
// // dof0[1] = RandomGenerator_MT19937::genrand_real3();
// // dof0[2] = RandomGenerator_MT19937::genrand_real3();
//       nb_pair_force <3, double,
//       		     nb_interaction_geometric_none_tag,
//       		     nb_interaction_accelteration_none_tag,
//       		     nb_interaction_electrostatic_null_tag,
//       		     nb_interaction_vanderWaals_cutoff_tag,
//       		     nb_interaction_compute_f_tag>
//       	  (eleParam, vdwParam, cutoff, &index0, &index1, dof, charge,
//       	   force, NULL, NULL);
//       // double diff[MDDIM];
//       // nb_auxiliary_diff <MDDIM, double> (dof0, dof1, diff);
//       // double dist2 = nb_auxiliary_dist2<MDDIM, double>  (diff);
//       // if (dist2 > cutoff[0] * cutoff[0]) continue;
//       // double ri2 = 1./dist2;
//       // double sri2 = vdwParam[LennardJones6_12::sigma] * vdwParam[LennardJones6_12::sigma] * ri2;
//       // double sri6 = sri2*sri2*sri2;
//       // double scalor = 24 * vdwParam[LennardJones6_12::epsilon] * (double(2) * (sri6*sri6) - sri6) * ri2;      
//       // force[0] = scalor * diff[0];
//       // force[1] = scalor * diff[1];
//       // force[2] = scalor * diff[2];
//     }
//   }
//   // nb_pair_force <3, double,
//   // 		 nb_interaction_geometric_none_tag,
//   // 		 nb_interaction_accelteration_none_tag,
//   // 		 nb_interaction_electrostatic_null_tag,
//   // 		 nb_interaction_vanderWaals_cutoff_tag,
//   // 		 nb_interaction_compute_f_tag>
//   //     (eleParam, vdwParam, cutoff, NULL, NULL, dof0, dof1, &charge0, &charge1,
//   //      NULL, force, NULL);

//   cout << force[0] << endl;
//   cout << force[1] << endl;
//   cout << force[2] << endl;

//   return 0;
// }
