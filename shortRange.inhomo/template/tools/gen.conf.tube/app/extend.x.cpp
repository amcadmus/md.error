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
#include "RandomGenerator.h"
#include "GenTubeAux.h"
#include "GroFileManager.h"

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  double Lx;
  std::string filename_in;
  std::string filename_out;
  po::options_description desc ("Allow options");  

  desc.add_options()
      ("help,h", "print this message")
      ("Lx", po::value<double > (&Lx)->default_value (150.),  "box size x")
      ("file-name-in",  po::value<std::string > (&filename_in)-> default_value (std::string("conf.gro"),  "input file"))
      ("file-name-out",	po::value<std::string > (&filename_out)->default_value (std::string("confout.gro"), "output file"));
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::vector<double > box(3);
  std::vector<std::vector<double > > posi, velo;
  std::vector<double > zero (3, 0.);
  std::vector<std::string > resdname, atomname;
  std::vector<int > resdindex, atomindex;
  GroFileManager::read (filename_in,
			resdindex, resdname,
			atomname, atomindex,
			posi, velo, box);
  if (box[0] > Lx){
    std::cerr << "new box Lx should be larger than current box Lx" << std::endl;
    return 1;
  }
  
  double shift = (Lx - box[0]) * 0.5;
  std::cout << "# old Lx is " << box[0] << std::endl;
  std::cout << "# new Lx is " << Lx << std::endl;
  std::cout << "# shift is " << shift << std::endl;
  std::cout << "# out file is " << filename_out << std::endl;

  for (unsigned i = 0; i < posi.size(); ++i){
    if      (posi[i][0] <  box[0]) posi[i][0] += box[0];
    else if (posi[i][0] >= box[0]) posi[i][0] -= box[0];
    posi[i][0] += shift;
  }

  box[0] = Lx;
  GroFileManager::write (filename_out,
			 resdindex, resdname,
			 atomname, atomindex,
			 posi, velo, box);
  return 0;
}

