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
#include "GroFileManager.h"

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  double L, rball, refh;
  
  std::string filename;
  po::options_description desc ("Allow options");  

  desc.add_options()
      ("help,h", "print this message")
      ("box-size,L", po::value<double > (&L)    ->default_value (84.553),  "box size.")
      ("lattice,h",  po::value<double > (&refh) ->default_value (1.5658),  "lattice size.")
      ("r-ball,r",   po::value<double > (&rball)->default_value (13.5)  ,  "radius of the cluster.")
      ("file-name,f",	po::value<std::string > (&filename)->	default_value (std::string("confout.gro"), "file to output"));
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::vector<std::vector<double > > posi, velo;
  std::vector<double > zero (3, 0.);
  std::vector<std::string > resdname, atomname;
  std::vector<int > resdindex, atomindex;
  std::vector<double > box(3);
  box[0] = L;
  box[1] = L;
  box[2] = L;

  int nlattice = int(L / refh) + 1;
  refh = L / (nlattice);
  std::cout << "# lattice size is " << refh << std::endl;
  std::cout << "# nlattice is     " << nlattice << std::endl;
  std::cout << "# box is          " << L << std::endl;
  
  for (int ix = 0; ix < nlattice; ++ix){
    for (int iy = 0; iy < nlattice; ++iy){
      for (int iz = 0; iz < nlattice; ++iz){
	std::vector<double > tmp(3);
	tmp[0] = (ix + 0.25) * refh;
	tmp[1] = (iy + 0.25) * refh;
	tmp[2] = (iz + 0.25) * refh;
	for (unsigned i = 0; i < 4; ++i){
	  if (i != 3){
	    if (i != 0) tmp[0] += 0.5 * refh;
	    if (i != 1) tmp[1] += 0.5 * refh;
	    if (i != 2) tmp[2] += 0.5 * refh;
	  }
	  double drx = tmp[0] - L * 0.5;
	  double dry = tmp[1] - L * 0.5;
	  double drz = tmp[2] - L * 0.5;
	  double diff = sqrt (drx * drx + dry * dry + drz * drz);
	  if (diff < rball){
	    posi.push_back (tmp);
	  }
	}
      }
    }
  }

  unsigned natom = posi.size();
  std::cout << "# natom is        " << natom << std::endl;
  unsigned count = 1;
  for (unsigned i = 0; i < natom; ++i){
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

