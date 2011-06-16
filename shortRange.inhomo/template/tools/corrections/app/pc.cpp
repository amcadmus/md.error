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
  std::string filename;
  po::options_description desc ("Allow options");  

  float start_t, end_t;
  double h;
  double rcut;
  
  desc.add_options()
      ("help,h", "print this message")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("rcut,r",  po::value<double > (&rcut)->default_value (5.),  "rcut-off")
      ("bin-size",  po::value<double > (&h)->default_value (1.),  "bin size")
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("traj.xtc"), "trajactory file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  DensityProfile_PiecewiseConst dp (filename, h);
  dp.deposit (filename, start_t, end_t);
  dp.print_x ("density.x.out");

  PresureCorrection pc (rcut, dp);
  pc.correction (dp);
  printf ("pc is %e %e %e\n", pc.pxx, pc.pyy, pc.pzz);
  printf ("tension is %e\n", (pc.pxx - (pc.pyy + pc.pzz) * 0.5) * 0.5 * dp.getBox()[0]);
  
  return 0;
}

 
 
