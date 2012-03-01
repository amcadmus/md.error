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
  std::vector<double > shift (3);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("shift-x,x", po::value<double > (&shift[0])->default_value(0.), "shift on x")
      ("shift-y,y", po::value<double > (&shift[1])->default_value(0.), "shift on y")
      ("shift-z,z", po::value<double > (&shift[2])->default_value(0.), "shift on z")
      ("input,i",   po::value<std::string > (&ifile)->default_value ("conf.gro"), "the input file")
      ("output,o",  po::value<std::string > (&ofile)->default_value ("confout.gro"), "the output file");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  printf ("#######################################################\n");
  printf ("## shift: %f %f %f\n", shift[0], shift[1], shift[2]);
  printf ("#######################################################\n");

  std::vector<int > resdindex, atomindex;
  std::vector<std::string > resdname, atomname;
  std::vector<std::vector<double > > posi, velo;
  std::vector<double > box;
  
  GroFileManager::read (ifile, resdindex, resdname, atomname, atomindex,
			posi, velo, box);
  for (unsigned i = 0; i < posi.size(); ++i){
    for (unsigned dd = 0; dd < 3; ++dd){
      posi[i][dd] += shift[dd];
      if (posi[i][dd] >= box[dd]){
	posi[i][dd] -= box[dd];
      }
      else if (posi[i][dd] < 0.){
	posi[i][dd] += box[dd];
      }
    }
  }
  GroFileManager::write (ofile, resdindex, resdname, atomname, atomindex,
			 posi, velo, box);
  
  return 0;
}

