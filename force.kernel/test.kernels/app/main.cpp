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

#include "nb_kernels.h"

int main(int argc, char * argv[])
{
  double eleParam[5] = {0.};
  double vdwParam[5] = {0.};
  double cutoff[1] = {1.0};
  
  
  
  nb_pair_force <3, double,
		 nb_interaction_geometric_none_tag,
		 nb_interaction_accelteration_none_tag,
		 nb_interaction_electrostatic_null_tag,
		 nb_interaction_vanderWaals_cutoff_tag,
		 nb_interaction_compute_f_tag>
      (eleParam, vdwParam, cutoff, NULL, NULL, NULL, NULL, NULL, NULL,
       NULL, NULL, NULL);
  
}
