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
#include <time.h>

#include <boost/program_options.hpp>
#include "DensityProfile.h"
#include "PresureCorrection.h"

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  std::string density_save;
  std::string rcut_save;
  po::options_description desc ("Allow options");
  
  desc.add_options()
      ("help,h", "print this message")
      ("density-save-name,r", po::value<std::string > (&density_save)->default_value (std::string("density.save"), "density save file"))
      ("rcut-save-name,d",    po::value<std::string > (&rcut_save   )->default_value (std::string("rcut.save"   ), "rcut save file"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  DensityProfile_PiecewiseConst dp;
  dp.load (density_save.c_str());
  dp.print_x ("density.x.out");
  RCutTable rctable;
  rctable.load_rc (rcut_save.c_str());

  PresureCorrection_RCutTable pc (rctable, dp);
  pc.correction (rctable, dp);
  printf ("pc is %e %e %e\n", pc.pxx, pc.pyy, pc.pzz);
  printf ("tension is %e\n", (pc.pxx - (pc.pyy + pc.pzz) * 0.5) * 0.5 * dp.getBox()[0]);

  for (unsigned i = 0; i < rctable.getProfileIndex().size(); ++i){
    rctable.getProfileIndex()[i] = 9;
  }
  pc.correction (rctable, dp);
  printf ("pc is %e %e %e\n", pc.pxx, pc.pyy, pc.pzz);
  printf ("tension is %e\n", (pc.pxx - (pc.pyy + pc.pzz) * 0.5) * 0.5 * dp.getBox()[0]);
  
  return 0;
}

 
 
