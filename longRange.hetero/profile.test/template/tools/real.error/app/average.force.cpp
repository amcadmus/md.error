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

#include "GroFileManager.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  po::options_description desc ("Allow options");
  std::string force0, force1, out, config;
  
  desc.add_options()
      ("help,h", "print this message")
      ("config,c", po::value<std::string > (&config)->default_value("conf.gro"), "config file")
      ("force-0", po::value<std::string > (&force0)->default_value("force0.out"), "force file 0")
      ("force-1", po::value<std::string > (&force1)->default_value("force1.out"), "force file 1")
      ("output,o", po::value<std::string > (&out)->default_value("avg.force.out"), "error file")
      ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  FILE * fpout;
  fpout = fopen (out.c_str(), "w");
  if (fpout == NULL){
    std::cerr << "cannot open file " << out << std::endl;
    return 1;
  }

  FILE * fp0, * fp1;
  fp0 = fopen (force0.c_str(), "r");
  fp1 = fopen (force1.c_str(), "r");
  if (fp0 == NULL){
    std::cerr << "cannot open file " << force0 << std::endl;
    return 1;
  }
  if (fp1 == NULL){
    std::cerr << "cannot open file " << force1 << std::endl;
    return 1;
  }
  
  while (1){
    double dir_f0x, dir_f0y, dir_f0z;
    double dir_f1x, dir_f1y, dir_f1z;
    double rec_f0x, rec_f0y, rec_f0z;
    double rec_f1x, rec_f1y, rec_f1z;
    double t_f0x, t_f0y, t_f0z;
    double t_f1x, t_f1y, t_f1z;
    int count0, count1;
    count0 = fscanf (fp0, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     &dir_f0x, &dir_f0y, &dir_f0z,
		     &rec_f0x, &rec_f0y, &rec_f0z,
		     &t_f0x, &t_f0y, &t_f0z
	);
    count1 = fscanf (fp1, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     &dir_f1x, &dir_f1y, &dir_f1z,
		     &rec_f1x, &rec_f1y, &rec_f1z,
		     &t_f1x, &t_f1y, &t_f1z
	);
    if (count0 == EOF || count1 == EOF){
      break;
    }
    else if (count0 != 9 || count1 != 9){
      std::cerr << "format error of the force file" << std::endl;
      break;
    }

    fprintf (fpout, "%.16e %.16e %.16e   %.16e %.16e %.16e   %.16e %.16e %.16e\n",
	     0.5 * (dir_f0x + dir_f1x),
	     0.5 * (dir_f0y + dir_f1y),
	     0.5 * (dir_f0z + dir_f1z),
	     0.5 * (rec_f0x + rec_f1x),
	     0.5 * (rec_f0y + rec_f1y),
	     0.5 * (rec_f0z + rec_f1z),
	     0.5 * (t_f0x + t_f1x),
	     0.5 * (t_f0y + t_f1y),
	     0.5 * (t_f0z + t_f1z));
  }
  
  fclose (fp0);
  fclose (fp1);
  fclose (fpout);  
}


