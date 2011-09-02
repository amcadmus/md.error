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
#include "AssignFcorr.h"

namespace po = boost::program_options;

#include <stdio.h>

int main (int argc, char * argv[])
{
  std::string filename;
  std::string filea, fileb, filed;
  po::options_description desc ("Allow options");  

  float start_t, end_t;
  
  desc.add_options()
      ("help,h", "print this message")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("fc.ctj"), "rcut file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::cout << "# start at: " << start_t << std::endl;
  std::cout << "# end   at: " << end_t   << std::endl;

  AssignFcorr assign;
  float time;
  float time_prec = 0.001;
  assign.init_read (filename.c_str());
  while (assign.read(time)){
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
    char timefile[1024];
    int zhengshu = int(time);
    int xiaoshu  = int((time - zhengshu) * 100.);
    sprintf (timefile, "fcorr.t%04d.%02d.out", zhengshu, xiaoshu);
    assign.print_xy (timefile);
  }
  
  return 0;
}

