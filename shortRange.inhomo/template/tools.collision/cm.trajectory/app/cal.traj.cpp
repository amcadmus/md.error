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
#include "DensityProfile.h"

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  std::string filename;
  std::string filea, fileb, filed;
  po::options_description desc ("Allow options");  

  float start_t, end_t;
  double h;
  
  desc.add_options()
      ("help,h", "print this message")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("bin-size,b",  po::value<double > (&h)->default_value (1.),  "bin size")
      ("out-file-a",  po::value<std::string > (&filea)->default_value ("traj.a.out"),  "trajectory of center of mass cluster a")
      ("out-file-b",  po::value<std::string > (&fileb)->default_value ("traj.b.out"),  "trajectory of center of mass cluster b")
      ("out-file-dist",  po::value<std::string > (&filed)->default_value ("dist.cm.out"),  "distance of the centers of mass")
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("traj.xtc"), "trajactory file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::cout << "# start at: " << start_t << std::endl;
  std::cout << "# end   at: " << end_t   << std::endl;
  FILE * fpa = fopen (filea.c_str(), "w");
  FILE * fpb = fopen (fileb.c_str(), "w");
  FILE * fpd = fopen (filed.c_str(), "w");
  if (fpa == NULL){
    std::cerr << "cannot open file " << filea << std::endl;
    return 1;
  }
  if (fpb == NULL){
    std::cerr << "cannot open file " << fileb << std::endl;
    return 1;
  }
  if (fpd == NULL){
    std::cerr << "cannot open file " << filed << std::endl;
    return 1;
  }

  double cma[3] = {0,0,0};
  double cmb[3] = {0,0,0};
  XDRFILE * fp = xdrfile_open (filename.c_str(), "r");
  int natoms;
  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  read_xtc_natoms (fname, &natoms);
  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  float time_prec = 0.001;
  while (read_xtc (fp, natoms, &step, &time, gbox, xx, &prec) == 0){
    if (end_t != 0.f) {
      if (time < start_t - time_prec){
	continue;
      }
      else if (time > end_t + time_prec) {
	break;
      }	
    }
    else {
      if (time < start_t - time_prec) continue;
    }
    std::cout << "#! loaded frame at time " << time << "ps   \r";  
    std::cout << std::flush;
    //////////////////////////////////////////////////
    // calculate the cms
    cma[2] = cma[1] = cma[0] = 0.;
    cmb[2] = cmb[1] = cmb[0] = 0.;
    for (int i = 0; i < natoms/2; ++i){
      cma[0] += xx[i][0];
      cma[1] += xx[i][1];
      cma[2] += xx[i][2];
    }
    for (int i = natoms/2; i < natoms; ++i){
      cmb[0] += xx[i][0];
      cmb[1] += xx[i][1];
      cmb[2] += xx[i][2];
    }
    cma[0] /= double (natoms/2);
    cma[1] /= double (natoms/2);
    cma[2] /= double (natoms/2);
    cmb[0] /= double (natoms/2);
    cmb[1] /= double (natoms/2);
    cmb[2] /= double (natoms/2);
    double dx = cma[0] - cmb[0];
    double dy = cma[1] - cmb[1];
    double dz = cma[2] - cmb[2];
    double diff = sqrt (dx*dx + dy*dy + dz*dz);
    fprintf (fpa, "%f %f %f %f\n", cma[0], cma[1], cma[2], time);
    fprintf (fpb, "%f %f %f %f\n", cmb[0], cmb[1], cmb[2], time);
    fprintf (fpd, "%f %f\n", time, diff);
  }
  printf ("\n");

  fclose (fpa);
  fclose (fpb);
  fclose (fpd);
  free (xx);
  
  return 0;
}

  
