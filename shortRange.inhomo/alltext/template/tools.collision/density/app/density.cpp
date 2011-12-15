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

  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  XDRFILE * fp = xdrfile_open (fname, "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return 1;
  }
  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  int natoms = dp.getNatoms();
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  float time_prec = 0.001;

  while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
    if (end_t != 0.f) {
      if (time < start_t - time_prec){
	continue;
      }
      else if (time > end_t + time_prec) {
	break;
      }	
    }
    else {
      if (time < start_t - time_prec) continue;
    }
    std::cout << "#! loaded frame at time " << time << "ps   \r";  
    std::cout << std::flush;
    dp.deposit (xx, natoms);
    char timefile[1024];
    int zhengshu = int(time);
    int xiaoshu  = int((time - zhengshu) * 100.);
    sprintf (timefile, "rho.t%04d.%02d.out", zhengshu, xiaoshu);
    dp.print_xy (timefile);
    dp.clear();
  }
  
  std::cout << std::endl;

  free(xx);
  xdrfile_close(fp);

  return 0;
}

  
