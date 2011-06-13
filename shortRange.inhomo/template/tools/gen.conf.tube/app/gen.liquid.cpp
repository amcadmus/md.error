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
  double Ly, Lz;
  double rho;
  int natom;
  
  std::string filename;
  po::options_description desc ("Allow options");  

  desc.add_options()
      ("help,h", "print this message")
      ("natom,n", po::value<int > (&natom)->default_value (16000),  "number of atom")
      ("Ly", po::value<double > (&Ly)->default_value (20.),  "box size y")
      ("Lz", po::value<double > (&Lz)->default_value (20.),  "box size z")
      ("rho", po::value<double > (&rho)->default_value (0.85),  "higher density")
      ("file-name,f",	po::value<std::string > (&filename)->	default_value (std::string("confout.gro"), "file to output"));
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  Lx = natom / rho / Ly / Lz;
  
  unsigned seed = unsigned (time(NULL));
  std::cout << "# seed is " << seed << std::endl;
  std::cout << "# Lx is " << Lx << std::endl;
  std::cout << "# Ly is " << Ly << std::endl;
  std::cout << "# Lz is " << Lz << std::endl;
  std::cout << "# rho is " << rho << std::endl;
  std::cout << "# out file is " << filename << std::endl;
  
  RandomGenerator_MT19937::init_genrand(seed);

  std::vector<std::vector<double > > posi, velo;
  std::vector<double > zero (3, 0.);
  std::vector<std::string > resdname, atomname;
  std::vector<int > resdindex, atomindex;
  unsigned count = 1;
  std::vector<double > box(3);
  box[0] = Lx;
  box[1] = Ly;
  box[2] = Lz;
  
  for (int i = 0; i < natom; ++i){
    std::vector<double > tmp(3);
    tmp[0] = Lx * RandomGenerator_MT19937::genrand_real2();
    tmp[1] = Ly * RandomGenerator_MT19937::genrand_real2();
    tmp[2] = Lz * RandomGenerator_MT19937::genrand_real2();
    
    posi.push_back(tmp);
    velo.push_back(zero);
    resdname.push_back(std::string("lj"));
    atomname.push_back(std::string("lj"));
    resdindex.push_back(count);
    atomindex.push_back(count);
    count ++;
  }
  
  GroFileManager::write (filename,
			 resdindex, resdname,
			 atomname, atomindex,
			 posi, velo, box);
  return 0;
}

