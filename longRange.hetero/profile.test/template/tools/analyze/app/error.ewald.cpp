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
#include "ErrorEstimate_Ewald.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  std::string tfile, efile, rfile, qfile;
  double beta;
  IntVectorType Kmax;
  IntVectorType K;
  int kValue, kmaxValue;
  float start, end;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("trajactory,t", po::value<std::string > (&tfile)->default_value("traj.xtc"), "trajactory file")
      ("charge-table,q", po::value<std::string > (&qfile), "charge table")      
      ("start,s", po::value<float > (&start)->default_value(0.f), "start time")
      ("end,e",   po::value<float > (&end  )->default_value(0.f), "end   time")
      ("beta,b", po::value<double > (&beta)->default_value(1.), "value of beta")
      ("kx", po::value<int > (&K.x)->default_value (27), "Number of grid points, should be odd")
      ("ky", po::value<int > (&K.y)->default_value (27), "Number of grid points, should be odd")
      ("kz", po::value<int > (&K.z)->default_value (27), "Number of grid points, should be odd")
      ("grid,k",po::value<int > (&kValue), "Number of grid points, should be odd, this will overwrite kx, ky and kz")
      ("kmaxx", po::value<int > (&Kmax.x)->default_value (9), "Cutoff in reciprocal space, should be odd")
      ("kmaxy", po::value<int > (&Kmax.y)->default_value (9), "Cutoff in reciprocal space, should be odd")
      ("kmaxz", po::value<int > (&Kmax.z)->default_value (9), "Cutoff in reciprocal space, should be odd")
      ("kmax",  po::value<int > (&kmaxValue), "Cutoff in reciprocal space, should be odd, this will overwrite kmaxx, kmaxy and kmaxz")
      ("output-density",  po::value<std::string > (&rfile)->default_value ("rho.x.avg.out"), "the output density (averaged on yz) of the system")
      ("output-error,o",  po::value<std::string > (&efile)->default_value ("error.out"), "the output error of the system");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  if (vm.count("grid")){
    K.z = K.y = K.x = kValue;
  }
  if (vm.count("kmax")){
    Kmax.z = Kmax.y = Kmax.x = kmaxValue;
  }

  printf ("#######################################################\n");
  printf ("## start traj at %f\n", start);
  printf ("## end   traj at %f\n", end);
  printf ("## beta is %.2f\n", beta);  
  printf ("## K    is %d %d %d\n", K.x, K.y, K.z);
  printf ("## Kmax is %d %d %d\n", Kmax.x, Kmax.y, Kmax.z);
  if (vm.count("charge-table")){
    printf ("## charge table: %s", qfile.c_str());
  }
  printf ("#######################################################\n");
  

  DensityProfile_PiecewiseConst dp;

  if (vm.count("charge-table")){
    dp.reinit_xtc (tfile.c_str(), qfile.c_str(), K.x, K.y, K.z, start, end);
  }
  else {
    dp.reinit_xtc (tfile.c_str(), K.x, K.y, K.z, start, end);
  }
  dp.print_avg_x (rfile.c_str());

  ErrorEstimate_Ewald eee;
  eee.reinit (beta, Kmax, dp);
  eee.calError (dp);
  eee.print_error (efile.c_str());
  eee.print_meanf ("meanf.out", dp);

  return 0;
}
