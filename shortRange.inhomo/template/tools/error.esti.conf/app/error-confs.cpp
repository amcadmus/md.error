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

#include "ErrorEstimator.h"
#include "DensityProfile.h"

int main(int argc, char * argv[])
{
  if (argc != 5){
    std::cerr << "usage:\n"
	     << argv[0] << " fileFilenames refh_Profile refh_Error rcut" << std::endl;
    return 1;
  }
  double refh_Profile, refh_Error, rcut;
  refh_Profile = atof(argv[2]);
  refh_Error = atof(argv[3]);
  rcut = atof(argv[4]);
  std::vector<std::string> names;
  FILE * fp = fopen (argv[1], "r");
  char name[1024];
  int c;
  while ((c = fscanf (fp, "%s", name)) == 1) {
    names.push_back (std::string(name));
    std::cout << "# find file " << std::string(name) << " c is " << c << std::endl;
  }
  fclose (fp);

  DensityProfile_PiecewiseConst dp (names, refh_Profile);
  dp.print_x  (std::string("profile.x.out"));
  dp.print_xy (std::string("profile.xy.out"));
  
  Disperson6 lj (1., 1., rcut);
  ErrorEstimatorConvN2 error (lj, dp.getBox(), refh_Error);

  error.estimate (dp);
  error.print_x (std::string("esti.x.out"));
  error.print_xy(std::string("esti.xy.out"));
  
  return 0;
}


