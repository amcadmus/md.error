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
#include "AdaptRCut.h"

int main(int argc, char * argv[])
{  
  if (argc != 3){
    std::cerr << "usage:\n"
	      << argv[0] << " xtc refh_Profile" << std::endl;
    return 1;
  }
  double refh_Profile;
  refh_Profile = atof(argv[2]);

  DensityProfile_Corr_PiecewiseConst dp (std::string(argv[1]), 0, refh_Profile);
  dp.print_x  (std::string("profile.x.out"));
  dp.print_xy (std::string("profile.xy.out"));

  AdaptRCut adprcut;
  adprcut.reinit (3., 10., 0.5, dp);
  adprcut.calError (dp);
  adprcut.print_x  (std::string("esti.x.out"));
  adprcut.calRCut  (2e-3);
  adprcut.print_rc (std::string("rcut.x.out"));
  

  // error.estimate (dp, rcorr * ratio, rcut);
  // error.print_x (std::string("esti.x.out"));
  // error.print_xy(std::string("esti.xy.out"));
  
  return 0;
}


