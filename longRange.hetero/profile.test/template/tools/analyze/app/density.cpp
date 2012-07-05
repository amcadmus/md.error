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

#include "DensityProfile.h"
#include "ErrorEstimate_dir.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void decideK (const std::string & fname,
	      const double & refh,
	      IntVectorType & K)
{
  XDRFILE * fp = xdrfile_open (fname.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << fname << std::endl;
    exit(1);
  }
  char tfname [1024];
  int natoms;
  strncpy (tfname, fname.c_str(), 1024);
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
  if (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) != 0){
    std::cerr << "error while loading xtc file" << std::endl;
    exit(1); 
  }

  std::vector<double > boxsize(3);
  boxsize[0] = gbox[0][0];
  boxsize[1] = gbox[1][1];
  boxsize[2] = gbox[2][2];
  K.x = unsigned (boxsize[0] / refh);
  K.y = unsigned (boxsize[1] / refh);
  K.z = unsigned (boxsize[2] / refh);
  if ((K.x - (K.x/2)*2) == 0) K.x++;
  if ((K.y - (K.y/2)*2) == 0) K.y++;
  if ((K.z - (K.z/2)*2) == 0) K.z++;

  free (xx);
  xdrfile_close (fp);
}


int main(int argc, char * argv[])
{
  std::string tfile, rfile, qfile;
  double refh;
  IntVectorType K;
  float start, end;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("start,s", po::value<float > (&start)->default_value(0.f), "start time")
      ("end,e",   po::value<float > (&end  )->default_value(0.f), "end   time")
      ("trajactory,t", po::value<std::string > (&tfile)->default_value("traj.xtc"), "trajactory file")
      ("charge-table,q", po::value<std::string > (&qfile)->default_value("charge.tab"), "charge table")
      ("refh,n", po::value<double > (&refh)->default_value(0.5), "bin size")
      ("output-density",  po::value<std::string > (&rfile)->default_value ("dir.avg.out"), "the output direction (averaged on yz) of the system");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  decideK (tfile, refh, K);

  printf ("#######################################################\n");
  printf ("## start traj at %f\n", start);
  printf ("## end   traj at %f\n", end);
  printf ("## K    is %d %d %d\n", K.x, K.y, K.z);
  if (vm.count("charge-table")){
    printf ("## charge table: %s", qfile.c_str());
  }
  printf ("#######################################################\n");

  DensityProfile_PiecewiseConst dp;

  dp.reinit_xtc (tfile.c_str(), qfile.c_str(), K.x, K.y, K.z, start, end);
  dp.print_avg_x (rfile.c_str());

  return 0;
}
