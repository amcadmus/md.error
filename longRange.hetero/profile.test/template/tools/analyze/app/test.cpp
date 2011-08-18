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

#include "DensityProfile.h"
#include "ErrorEstimate_Ewald.h"

int main(int argc, char * argv[])
{
  DensityProfile_PiecewiseConst dp;

  // dp.reinit_conf ("conf.gro", 81, 81, 81);
  dp.reinit_xtc ("traj.xtc", 81, 81, 81, 0, 0);
  dp.print_avg_x ("rho.x.avg.out");

  ErrorEstimate_Ewald eee;
  IntVectorType Kmax;
  Kmax.x = Kmax.y = Kmax.z = 9;
  eee.reinit (1.0, Kmax, dp);
  eee.calError (dp);
  eee.print_error ("error.out");
}
