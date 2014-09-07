#include "DensityFunction.h"
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
#include "RandomGenerator.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "GroFileManager.h"

int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  std::vector<double > cut (3);
  double extension = 0.;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("new-x,x", po::value<double > (&cut[0])->default_value(-1.), "cut on x")
      ("new-y,y", po::value<double > (&cut[1])->default_value(-1.), "cut on y")
      ("new-z,z", po::value<double > (&cut[2])->default_value(-1.), "cut on z")
      ("extension,e", po::value<double > (&extension)->default_value(0.1), "extension")
      ("input,i",   po::value<std::string > (&ifile)->default_value ("conf.gro"), "the input file")
      ("output,o",  po::value<std::string > (&ofile)->default_value ("out.gro"), "the output file");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  printf ("#######################################################\n");
  printf ("## cut: %f %f %f\n", cut[0], cut[1], cut[2]);
  printf ("#######################################################\n");

  std::vector<int > resdindex, atomindex;
  std::vector<std::string > resdname, atomname;
  std::vector<std::vector<double > > posi, velo;
  std::vector<double > box;
  
  GroFileManager::read (ifile, resdindex, resdname, atomname, atomindex,
			posi, velo, box);

  std::vector<int > resdindex1, atomindex1;
  std::vector<std::string > resdname1, atomname1;
  std::vector<std::vector<double > > posi1, velo1;
  std::vector<double > box1;

  for (unsigned ii = 0; ii < posi.size() / 3; ++ii){
    bool include = true;
    for (unsigned dd = 0; dd < 3; ++dd){
      if (cut[dd] > 0 && posi[ii*3][dd] >= cut[dd]) {
	include = false;
      }
    }
    if (include){
      for (unsigned jj = 0; jj < 3; ++jj){
	posi1.push_back (posi[ii*3+jj]);
	velo1.push_back (velo[ii*3+jj]);
	resdindex1.push_back (resdindex[ii*3+jj]);
	atomindex1.push_back (atomindex[ii*3+jj]);
	resdname1.push_back (resdname[ii*3+jj]);
	atomname1.push_back (atomname[ii*3+jj]);
      }
    }
  }
  
  box1 = box;
  for (unsigned dd = 0; dd < 3; ++dd){
    if (cut[dd] > 0){
      box1[dd] = cut[dd] + extension;
    }
  }
	
  GroFileManager::write (ofile, resdindex1, resdname1, atomname1, atomindex1,
			 posi1, velo1, box1);
  
  return 0;
}

