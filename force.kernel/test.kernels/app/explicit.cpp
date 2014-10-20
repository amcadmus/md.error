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
#include "nb_interaction_lj.h"
#include "RandomGenerator.h"

int main(int argc, char * argv[])
{
  double vdwParam[LennardJones6_12::nparam];
  vdwParam[LennardJones6_12::epsilon] = 1.0;
  vdwParam[LennardJones6_12::sigma] = 1.0;

  double dof0[MDDIM], dof1[MDDIM];
  double force[MDDIM];
  force[0] = force[1] = force[2] = 0.;

  int n0 = 1000;
  int n1 = 100000;

  for (int ii = 0; ii < n0; ++ii){
    for (int jj = 0; jj < n1; ++jj){
      dof1[0] = RandomGenerator_MT19937::genrand_real3();
      dof1[1] = RandomGenerator_MT19937::genrand_real3();
      dof1[2] = RandomGenerator_MT19937::genrand_real3();
      dof0[0] = RandomGenerator_MT19937::genrand_real3();
      dof0[1] = RandomGenerator_MT19937::genrand_real3();
      dof0[2] = RandomGenerator_MT19937::genrand_real3();

      double diff[MDDIM];
      diff[0] = dof1[0] - dof0[0];
      diff[1] = dof1[1] - dof0[1];
      diff[2] = dof1[2] - dof0[2];
      double dist2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
      double ri2 = 1./dist2;
      double sri2 = vdwParam[LennardJones6_12::sigma] * vdwParam[LennardJones6_12::sigma] * ri2;
      double sri6 = sri2*sri2*sri2;
      double scalor = 24 * vdwParam[LennardJones6_12::epsilon] * (double(2) * (sri6*sri6) - sri6) * ri2;
      
      force[0] = scalor * diff[0];
      force[1] = scalor * diff[1];
      force[2] = scalor * diff[2];
    }
  }

  cout << force[0] << endl;
  cout << force[1] << endl;
  cout << force[2] << endl;

  return 0;
}
