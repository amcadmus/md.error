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
  std::string record, tfile;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("record-file,r", po::value<std::string > (&record)->default_value("conf.record"), "file of filenames to be analyzed")
      ("output-traj,o",  po::value<std::string > (&tfile)->default_value ("traj.xtc"), "the output traj file")
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  FILE * fp = fopen (record.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << record << std::endl;
    return 1;
  }
  XDRFILE * ofp = xdrfile_open (tfile.c_str(), "w");

  char confFileName[1024];
  int step = 0;
  while (fscanf (fp, "%s\n", confFileName) != EOF){
    std::vector<int > resdindex, atomindex;
    std::vector<std::string > resdname, atomname;
    std::vector<std::vector<double > > posi, velo;
    std::vector<double > tmpbox;
    GroFileManager::read (confFileName,
			  resdindex, resdname, atomname, atomindex,
			  posi, velo,
			  tmpbox);
    
    int natoms = posi.size();
    matrix box;
    box[0][0] = box[0][1] = box[0][2] = 0.;
    box[1][0] = box[1][1] = box[1][2] = 0.;
    box[2][0] = box[2][1] = box[2][2] = 0.;
    box[0][0] = tmpbox[0];
    box[1][1] = tmpbox[1];
    box[2][2] = tmpbox[2];
    rvec * x = (rvec *) malloc (sizeof(rvec) * natoms);
    for (int i = 0; i < natoms; ++i){
      x[i][0] = posi[i][0];
      x[i][1] = posi[i][1];
      x[i][2] = posi[i][2];
    }
    write_xtc (ofp, natoms, step, step*1., box, x, 1000.);
    
    free (x);
    step ++;
  }
    
  xdrfile_close (ofp);
  fclose (fp);
  
  return 0;
}

