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
  std::string tfile, efile, rfile;
  double beta, rcut, refh;
  IntVectorType K;
  float start, end;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("trajactory,t", po::value<std::string > (&tfile)->default_value("traj.xtc"), "trajactory file")
      ("start,s", po::value<float > (&start)->default_value(0.f), "start time")
      ("end,e",   po::value<float > (&end  )->default_value(0.f), "end   time")
      ("beta,b", po::value<double > (&beta)->default_value(1.), "value of beta")
      ("rcut,r", po::value<double > (&rcut)->default_value(2.), "value of rcut")
      ("refh,n", po::value<double > (&refh)->default_value(1.), "bin size")
      ("output-density",  po::value<std::string > (&rfile)->default_value ("rho.x.avg.out"), "the output density (averaged on yz) of the system")
      ("output-error,o",  po::value<std::string > (&efile)->default_value ("error.out"), "the output error of the system");
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
  printf ("## beta is %.2f\n", beta);  
  printf ("## rcut is %.2f\n", rcut);  
  printf ("## K    is %d %d %d\n", K.x, K.y, K.z);
  printf ("#######################################################\n");

  DensityProfile_PiecewiseConst dp;

  dp.reinit_xtc (tfile.c_str(), K.x, K.y, K.z, start, end);
  dp.print_avg_x (rfile.c_str());

  ErrorEstimate_dir eee;
  eee.reinit (beta, rcut, dp);
  eee.calError (dp);
  eee.print_error (efile.c_str());
  eee.print_meanf ("meanf.out", dp);

  return 0;
}
