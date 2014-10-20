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
#include "RandomGenerator.h"

const int MDDIM = 3;

int main(int argc, char * argv[])
{
  double eleParam[5] = {0.};
  double vdwParam[LennardJones6_12::nparam];
  vdwParam[LennardJones6_12::epsilon] = 1.0;
  vdwParam[LennardJones6_12::sigma] = 1.0;
  
  double cutoff[1] = {5.0};
  double dof0[MDDIM], dof1[MDDIM];
  dof1[0] = 1.5;
  dof1[1] = 0.0;
  dof1[2] = 0.0;
  dof0[0] = 0.0;
  dof0[1] = 0.0;
  dof0[2] = 0.0;
  double charge0 = 1;
  double charge1 = -1;
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
      nb_pair_force <3, double,
		     nb_interaction_geometric_none_tag,
		     nb_interaction_accelteration_none_tag,
		     nb_interaction_electrostatic_null_tag,
		     nb_interaction_vanderWaals_cutoff_tag,
		     nb_interaction_compute_f_tag>
	  (eleParam, vdwParam, cutoff, NULL, NULL, dof0, dof1, &charge0, &charge1,
	   NULL, force, NULL);
    }
  }
  // nb_pair_force <3, double,
  // 		 nb_interaction_geometric_none_tag,
  // 		 nb_interaction_accelteration_none_tag,
  // 		 nb_interaction_electrostatic_null_tag,
  // 		 nb_interaction_vanderWaals_cutoff_tag,
  // 		 nb_interaction_compute_f_tag>
  //     (eleParam, vdwParam, cutoff, NULL, NULL, dof0, dof1, &charge0, &charge1,
  //      NULL, force, NULL);

  cout << force[0] << endl;
  cout << force[1] << endl;
  cout << force[2] << endl;

  return 0;
}
