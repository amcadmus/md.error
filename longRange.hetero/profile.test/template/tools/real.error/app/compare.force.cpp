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
  std::string config, force0, force1, out;
  double refh;
  
  desc.add_options()
      ("help,h", "print this message")
      ("config,c", po::value<std::string > (&config)->default_value("conf.gro"), "config file")
      ("refh", po::value<double > (&refh)->default_value(1.), "bin size")
      ("force-0", po::value<std::string > (&force0)->default_value("force0.out"), "force file 0")
      ("force-1", po::value<std::string > (&force1)->default_value("force1.out"), "force file 1")
      ("error-file,o", po::value<std::string > (&out)->default_value("error.out"), "error file")
      ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

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
  
  std::vector<int > resdindex, atomindex;
  std::vector<std::string > resdname, atomname;
  std::vector<std::vector<double > > posi, velo;
  std::vector<double > tmpbox;
  GroFileManager::read (config,
			resdindex, resdname, atomname, atomindex,
			posi, velo,
			tmpbox);
  
  int nbin = tmpbox[0] / refh + 0.5;
  double hx = tmpbox[0] / double(nbin);
  std::vector<int > count(nbin, 0);
  std::vector<double > value (nbin, 0.);
  for (unsigned i = 0; i < posi.size(); ++i){
    if (posi[i][0] >= tmpbox[0]) posi[i][0] -= tmpbox[0];
    else if (posi[i][0] < 0) posi[i][0] += tmpbox[0];
    int idx = posi[i][0] / hx;
    double f0x, f0y, f0z;
    double f1x, f1y, f1z;
    fscanf (fp0, "%lf%lf%lf", &f0x, &f0y, &f0z);
    fscanf (fp1, "%lf%lf%lf", &f1x, &f1y, &f1z);
    f0x -= f1x;
    f0y -= f1y;
    f0z -= f1z;
    value[idx] += f0x*f0x + f0y*f0y + f0z * f0z;
    count[idx] += 1;
  }
  for (int i = 0; i < nbin; ++i){
    if (count[i] != 0){
      value[i] = sqrt (value[i] / count[i]);
    }
  }
  
  fclose (fp0);
  fclose (fp1);
  
  FILE * fp = fopen(out.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << out << std::endl;
    return 1;
  }
  for (int i = 0; i < nbin; ++i){
    fprintf(fp, "%f %e\n",
	    (i+0.5) * refh,
	    value[i]);
  }
  fclose (fp);
}

