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
#include "ErrorEstimate_SPME_Ik.h"
#include "ErrorEstimate_SPME_Ana.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  std::string tfile, efile, rfile, qfile;
  double beta, pcharge;
  IntVectorType K;
  int order;
  int kValue;
  float start, end;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("trajactory,t", po::value<std::string > (&tfile)->default_value("traj.xtc"), "trajactory file")
      ("charge-table,q", po::value<std::string > (&qfile), "charge table")
      ("my-charge", po::value<double > (&pcharge)->default_value(1.), "point positive charge as testing charge")
      ("start,s", po::value<float > (&start)->default_value(0.f), "start time")
      ("end,e",   po::value<float > (&end  )->default_value(0.f), "end   time")
      ("beta,b", po::value<double > (&beta)->default_value(1.), "value of beta")
      ("order,n", po::value<int > (&order)->default_value(4), "order of B-spline")
      ("kx", po::value<int > (&K.x)->default_value (27), "Number of grid points, should be odd")
      ("ky", po::value<int > (&K.y)->default_value (27), "Number of grid points, should be odd")
      ("kz", po::value<int > (&K.z)->default_value (27), "Number of grid points, should be odd")
      ("grid,k",po::value<int > (&kValue), "Number of grid points, should be odd, this will overwrite kx, ky and kz")
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

  printf ("#######################################################\n");
  printf ("## start traj at %f\n", start);
  printf ("## end   traj at %f\n", end);
  printf ("## beta  is %.2f\n", beta);  
  printf ("## K     is %d %d %d\n", K.x, K.y, K.z);
  printf ("## order is %d\n", order);
  printf ("## test charge is %.2f\n", pcharge);  
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

  ErrorEstimate_SPME_Ana eesi;
  // IntVectorType refine;
  // refine.x = refine.y = refine.z= 1;
  // refine.x = 8;
  // refine.y = 4;
  // refine.z = 4;
  eesi.reinit (beta, order, dp);
  eesi.calError (dp, pcharge);
  eesi.print_meanf ("meanf.out", dp);
  eesi.print_error ("error.out");

  return 0;
}
