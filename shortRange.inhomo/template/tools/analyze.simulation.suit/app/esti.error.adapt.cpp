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
  std::string rcutfilename;
  po::options_description desc ("Allow options");  
  
  desc.add_options()
      ("help,h", "print this message")
      ("ruct-file-name,s", po::value<std::string > (&rcutfilename)->default_value (std::string("rcut.save"), "rcut file name"))      
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("confout.gro"), "conf file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  int inttmp;
  double doubletmp;
  FILE * fp = fopen (rcutfilename.c_str(), "r");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", rcutfilename.c_str());
    exit (1);
  }
  fscanf (fp, "%d", &inttmp);
  for (int i = 0; i < inttmp; ++i){
    fscanf (fp, "%lf", &doubletmp);
  }
  double b, c;
  int nx, ny, nz;
  fscanf (fp, "%lf %lf %lf", &doubletmp, &b, &c);
  fscanf (fp, "%d %d %d", &nx, &ny, &nz);
  
  DensityProfile_PiecewiseConst dp;
  dp.reinit_conf (filename, nx, ny, nz);
  dp.print_x ("density.x.out");

  AdaptRCut arc;
  arc.load_rc (rcutfilename, dp);
  arc.print_error_avg (dp, "a1.error.x.out");

  return 0;
}

 
 
