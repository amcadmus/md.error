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

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  std::string filename;
  po::options_description desc ("Allow options");  

  float start_t, end_t;
  double h;
  
  desc.add_options()
      ("help,h", "print this message")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("bin-size,b",  po::value<double > (&h)->default_value (1.),  "bin size")
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
  dp.print_avg_x ("a.density.x.out");
  
  return 0;
}

  
