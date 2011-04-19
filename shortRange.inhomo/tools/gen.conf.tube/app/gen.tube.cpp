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
  // Lx: box length
  // Lh: high density region
  // Lt: transition layer thickness
  double Lx, Lh, Lt;
  double Ly, Lz;
  // densities
  double rhoh, rhol;
  
  std::string filename;
  po::options_description desc ("Allow options");  

  desc.add_options()
      ("help,h", "print this message")
      ("Lx", po::value<double > (&Lx)->default_value (60.),  "box size x")
      ("Ly", po::value<double > (&Ly)->default_value (20.),  "box size y")
      ("Lz", po::value<double > (&Lz)->default_value (20.),  "box size z")
      ("Lh", po::value<double > (&Lh)->default_value (20.),  "high density range")
      ("Lt", po::value<double > (&Lt)->default_value (20.),  "transition layer thickness")
      ("rhoh", po::value<double > (&rhoh)->default_value (1.),  "higher density")
      ("rhol", po::value<double > (&rhol)->default_value (.5),  "lower  density")
      ("file-name,f",	po::value<std::string > (&filename)->	default_value (std::string("confout.gro"), "file to output"));
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  double volume = Lx * Ly * Lz;
  unsigned natom = unsigned(volume / rhoh);
  unsigned seed = unsigned (time(NULL));
  RandomGenerator_MT19937::init_genrand(seed);
  AcceptRatio ar (Lx, Lh, Lt, rhoh, rhol);

  std::vector<std::vector<double > > posi, velo;
  std::vector<double > zero (3, 0.);
  std::vector<std::string > resdname, atomname;
  std::vector<int > resdindex, atomindex;
  unsigned count = 1;
  std::vector<double > box(3);
  box[0] = Lx;
  box[1] = Ly;
  box[2] = Lz;
  
  for (unsigned i = 0; i < natom; ++i){
    std::vector<double > tmp(3);
    tmp[0] = Lx * RandomGenerator_MT19937::genrand_real2();
    tmp[1] = Ly * RandomGenerator_MT19937::genrand_real2();
    tmp[2] = Lz * RandomGenerator_MT19937::genrand_real2();
    double value = ar (tmp[0]);
    
    if (RandomGenerator_MT19937::genrand_real1() <= value){
      posi.push_back(tmp);
      velo.push_back(zero);
      resdname.push_back(std::string("lj"));
      atomname.push_back(std::string("lj"));
      resdindex.push_back(count);
      atomindex.push_back(count);
      count ++;
    }
  }
  
  GroFileManager::write (filename,
			 resdindex, resdname,
			 atomname, atomindex,
			 posi, velo, box);
  return 0;
}

