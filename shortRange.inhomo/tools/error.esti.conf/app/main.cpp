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
  if (argc != 4){
    std::cerr << "usage:\n"
	     << argv[0] << " inputConf refh rcut" << std::endl;
    return 1;
  }

  DensityProfile_PiecewiseConst dp (std::string(argv[1]), atof(argv[2]));
  Disperson6 lj (1., 1., atof(argv[3]));
  ErrorEstimatorConvN2 error (lj);

  error.estimate (dp);
  error.print_x (std::string("esti.x.out"));
  error.print_xy(std::string("esti.xy.out"));
  
  return 0;
}


