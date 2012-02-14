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

#include "GroFileManager.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

int main(int argc, char * argv[])
{
  std::string ifile;
  int every;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("every,e", po::value<int > (&every)->default_value (1), "skip every this number of frames")
      ("input-traj,i", po::value<std::string > (&ifile)->default_value("traj.xtc"), "file of filenames to be analyzed")
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  int natoms;
  XDRFILE * fp = xdrfile_open (ifile.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << ifile << std::endl;
    exit(1);
  }
  char tfname [10240];
  strncpy (tfname, ifile.c_str(), 1024);
  if (read_xtc_natoms (tfname, &natoms) != 0){
    std::cerr << "wrong reading natoms" << std::endl;    
    exit (1);
  }
  std::cout << "# natom is " << natoms << std::endl;

  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * natoms);

  std::vector<double > box(3);
  std::vector<std::vector<double > > posi(natoms);
  std::vector<std::vector<double > > velo(natoms);
  std::vector<std::string > resdname(natoms);
  std::vector<std::string > atomname(natoms);
  std::vector<int > resdindex(natoms);
  std::vector<int > atomindex(natoms);
  for (int i = 0; i < natoms; ++i){
    resdname [i] = atomname [i] = std::string ("lj");
    resdindex[i] = atomindex[i] = i+1;
    std::vector<double > tmp (3, 0.);
    velo[i] = tmp;
  }
  int nfile = 0;
  int count = 0;
  while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
    if ((count ++) % every != 0) continue;
    nfile ++;
    char filename[10240];
    sprintf (filename, "conf.%04d.gro", nfile);
    std::string ofile (filename);
    std::cout << "# load frame " << count
	      << " at time " << time
	      << " ps"
	      << ", and write to file " << ofile << std::endl;
    std::vector<double > tmp (3, 0.);
    for (int i = 0; i < natoms; ++i){
      tmp[0] = xx[i][0];
      tmp[1] = xx[i][1];
      tmp[2] = xx[i][2];
      posi[i] = tmp;
    }
    box[0] = gbox[0][0];
    box[1] = gbox[1][1];
    box[2] = gbox[2][2];
    GroFileManager::write (ofile,
			   resdindex, resdname, atomname, atomindex,
			   posi, velo, box);
  }  
  xdrfile_close (fp);
  
  return 0;
}

