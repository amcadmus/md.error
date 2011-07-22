#define CPLUSPLUS

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
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include "xdrfile/xdrfile_trr.h"

namespace po = boost::program_options;

void crossProd (double *c,
		const double *a,
		const float *b)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}


int main (int argc, char * argv[])
{
  std::string filename;
  std::string fileout;
  po::options_description desc ("Allow options");  

  float start_t, end_t;
  
  desc.add_options()
      ("help,h", "print this message")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("out-file,o", po::value<std::string > (&fileout)->default_value (std::string("angluar.out"), "output file name"))
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("traj.trr"), "trajactory file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::cout << "# start at: " << start_t << std::endl;
  std::cout << "# end   at: " << end_t   << std::endl;

  FILE * fpo = fopen(fileout.c_str(), "w");
       
  XDRFILE * fp = xdrfile_open (filename.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return 1;
  }
  int natoms;
  char fname [1024];
  strncpy (fname, filename.c_str(), 1024);
  read_trr_natoms (fname, &natoms);
  int step;
  float time, lambda;
  matrix gbox;  
  rvec * xx, *vv, *ff;
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  vv = (rvec *) malloc (sizeof(rvec) * natoms);
  // ff = (rvec *) malloc (sizeof(rvec) * natoms);
  ff = NULL;
  float time_prec = 0.001;

  fprintf (fpo, "#    1   2   3   4 5   6   7   8\n");
  fprintf (fpo, "# time L_x L_y L_z I w_x w_y w_z\n");
  while (read_trr (fp, natoms, &step, &time, &lambda, gbox, xx, vv, ff) == 0){
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
    std::cout << "#! loaded frame at time " << time << "ps, ";  
    //////////////////////////////////////////////////
    // calculate the cms
    double com[3];
    com[2] = com[1] = com[0] = 0.;
    for (int i = 0; i < natoms; ++i){
      com[0] += xx[i][0];
      com[1] += xx[i][1];
      com[2] += xx[i][2];
    }
    com[0] /= double(natoms);
    com[1] /= double(natoms);
    com[2] /= double(natoms);
    double dev[3];
    dev[0] = com[0] - gbox[0][0] * 0.5;
    dev[1] = com[1] - gbox[1][1] * 0.5;
    dev[2] = com[2] - gbox[2][2] * 0.5;
    std::cout << "deviate from box center: "
	      << dev[0] << " "
	      << dev[1] << " "
	      << dev[2] << "            \r";
    std::cout << std::flush;
    double center[3];
    center[0] = gbox[0][0] * 0.5;
    center[1] = gbox[1][1] * 0.5;
    center[2] = gbox[2][0] * 0.5;
    double sum_moment[3];
    sum_moment[2] = sum_moment[1] = sum_moment[0] = 0.;
    double sum_rr = 0.;
    for (int i = 0; i < natoms; ++i){
      double dr[3];
      dr[0] = xx[i][0] - center[0];
      dr[1] = xx[i][1] - center[1];
      dr[2] = xx[i][2] - center[2];
      double mymoment[3];
      crossProd (mymoment, dr, vv[i]);
      double myrr = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      sum_moment[0] += mymoment[0];
      sum_moment[1] += mymoment[1];
      sum_moment[2] += mymoment[2];
      sum_rr += myrr;
    }
    fprintf (fpo, "%f  %e %e %e  %e  %e %e %e\n",
	     time,
	     sum_moment[0], 
	     sum_moment[1], 
	     sum_moment[2],
	     sum_rr,
	     sum_moment[0] / sum_rr,
	     sum_moment[1] / sum_rr,
	     sum_moment[2] / sum_rr);
  }
  printf ("\n");

  fclose (fpo);
  free (xx);
  free (vv);
  xdrfile_close (fp);
  // free (ff);
  
  return 0;
}

  
