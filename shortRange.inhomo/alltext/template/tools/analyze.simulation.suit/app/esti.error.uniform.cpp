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
#include "PressureCorrection.h"

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  std::string filename;
  po::options_description desc ("Allow options");  

  double h;
  double rcut;
  
  desc.add_options()
      ("help,h", "print this message")
      ("rcut,r", po::value<double >(&rcut)->default_value(5.), "cut-off radius")
      ("bin-size,b",  po::value<double > (&h)->default_value (1.),  "bin size")
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("confout.gro"), "conf file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  DensityProfile_PiecewiseConst dp;
  dp.reinit_conf (filename, h);
  // dp.print_x ("density.x.out");

  AdaptRCut arc;
  arc.reinit (rcut, rcut, 1, dp);
  arc.calError (dp);
  arc.uniformRCut (rcut);
  arc.print_error_avg (dp, "a.error.conf.x.out");

  return 0;
}

 
 
