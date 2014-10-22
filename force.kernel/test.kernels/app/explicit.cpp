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
const int MDDIM = 3;
typedef float ValueType;

#include "nb_interaction_lj.h"
#include "RandomGenerator.h"
#include "Stopwatch.h"

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
  // ValueType force0[MDDIM], force1[MDDIM];
  ValueType total[MDDIM] = {0.};

  Stopwatch mywatch;

  mywatch.start();
  for (int iindex = 0; iindex < n0; ++iindex) {
    ValueType dof0[MDDIM], dof1[MDDIM];
    dof0[0] = coord[iindex*3]	-2;
    dof0[1] = coord[iindex*3+1]	-2;
    dof0[2] = coord[iindex*3+2]	-2;
    ValueType tmpf [MDDIM] = {0.};
    
    for (int jj = neighbor_index_data[iindex*2]; jj < neighbor_index_data[iindex*2+1]; ++jj){
      int jindex = neighbor_index[jj];
      dof1[0] = coord[jindex*3];
      dof1[1] = coord[jindex*3+1];
      dof1[2] = coord[jindex*3+2];
      ValueType diff[MDDIM];
      diff[0] = dof0[0] - dof1[0];
      diff[1] = dof0[1] - dof1[1];
      diff[2] = dof0[2] - dof1[2];
      ValueType dist2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
      // if (dist2 > cutoff[0] * cutoff[0]) continue;
      ValueType ri2 = 1./dist2;
      ValueType sri2 = ri2;
      ValueType sri6 = sri2*sri2*sri2;
      ValueType scalor = (param[1] * (sri6*sri6) - param[0] * sri6) * ri2;
      
      // force0[0] = scalor * diff[0];
      // force0[1] = scalor * diff[1];
      // force0[2] = scalor * diff[2];
      force[jindex*3+0] +=-scalor * diff[0];
      force[jindex*3+1] +=-scalor * diff[1];
      force[jindex*3+2] +=-scalor * diff[2];
      tmpf[0] += scalor * diff[0];
      tmpf[1] += scalor * diff[1];
      tmpf[2] += scalor * diff[2];
    
      // total[0] += force1[0];
      // total[1] += force1[0];
      // total[2] += force1[0];
      // printf ("%f %f %f   %f\n", scalor * diff[0], scalor*diff[1], scalor*diff[2], scalor);
    }
    force[iindex*3+0] += tmpf[0];
    force[iindex*3+1] += tmpf[1];
    force[iindex*3+2] += tmpf[2];
  }
  
  mywatch.stop();

  for (int ii = 0; ii < n0; ++ii){
    printf ("%f %f %f\n", force[0+3*ii], force[1+3*ii], force[2+3*ii]);
  }

  printf ("syste: %f  user: %f  real: %f\n", mywatch.system(), mywatch.user(), mywatch.real());
  return 0;
}



// int main(int argc, char * argv[])
// {
//   ValueType vdwParam[LennardJones6_12::nparam];
//   vdwParam[LennardJones6_12::epsilon] = 1.0;
//   vdwParam[LennardJones6_12::sigma] = 1.0;
//   ValueType cutoff[1] = {5.0};

//   ValueType dof0[MDDIM], dof1[MDDIM];
//   dof1[0] = 1.5;
//   dof1[1] = 0.0;
//   dof1[2] = 0.0;
//   dof0[0] = 0.0;
//   dof0[1] = 0.0;
//   dof0[2] = 0.0;
//   ValueType force0[MDDIM];
//   force0[0] = force0[1] = force0[2] = 0.;
//   ValueType force1[MDDIM];
//   force1[0] = force1[1] = force1[2] = 0.;

//   int n0 = 1000;
//   int n1 = 100000;

//   for (int ii = 0; ii < n0; ++ii){
//     for (int jj = 0; jj < n1; ++jj){
//       dof1[0] = RandomGenerator_MT19937::genrand_real3();
//       // dof1[1] = RandomGenerator_MT19937::genrand_real3();
//       // dof1[2] = RandomGenerator_MT19937::genrand_real3();
//       // dof0[0] = RandomGenerator_MT19937::genrand_real3();
//       // dof0[1] = RandomGenerator_MT19937::genrand_real3();
//       // dof0[2] = RandomGenerator_MT19937::genrand_real3();

//       ValueType diff[MDDIM];
//       diff[0] = dof1[0] - dof0[0];
//       diff[1] = dof1[1] - dof0[1];
//       diff[2] = dof1[2] - dof0[2];
//       ValueType dist2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
//       if (dist2 > cutoff[0] * cutoff[0]) continue;
//       ValueType ri2 = 1./dist2;
//       ValueType sri2 = vdwParam[LennardJones6_12::sigma] * vdwParam[LennardJones6_12::sigma] * ri2;
//       ValueType sri6 = sri2*sri2*sri2;
//       ValueType scalor = 24 * vdwParam[LennardJones6_12::epsilon] * (ValueType(2) * (sri6*sri6) - sri6) * ri2;
      
//       force0[0] = scalor * diff[0];
//       force0[1] = scalor * diff[1];
//       force0[2] = scalor * diff[2];
//       force1[0] =-scalor * diff[0];
//       force1[1] =-scalor * diff[1];
//       force1[2] =-scalor * diff[2];
//     }
//   }

//   cout << force0[0] << endl;
//   cout << force0[1] << endl;
//   cout << force0[2] << endl;

//   return 0;
// }
