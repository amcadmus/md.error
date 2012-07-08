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
#include "GroFileManager.h"
#include "ErrorEstimate_SPME_St_H2O_gr.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  std::string tfile, efile, rfile, qfile, ofile;
  double beta, pcharge;
  IntVectorType K;
  int order;
  int kValue;
  float start, end;
  std::string fileOO, fileOH, fileHH;
  double rhoMol, grup;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("trajactory,t", po::value<std::string > (&tfile)->default_value("traj.xtc"), "trajactory file")
      ("charge-table,q", po::value<std::string > (&qfile)->default_value("charge.tab"), "charge table")
      ("my-charge", po::value<double > (&pcharge)->default_value(1.), "point positive charge as testing charge")
      ("start,s", po::value<float > (&start)->default_value(0.f), "start time")
      ("end,e",   po::value<float > (&end  )->default_value(0.f), "end   time")
      ("beta,b", po::value<double > (&beta)->default_value(1.), "value of beta")
      ("order,n", po::value<int > (&order)->default_value(4), "order of B-spline")
      ("kx", po::value<int > (&K.x)->default_value (27), "Number of grid points, should be odd")
      ("ky", po::value<int > (&K.y)->default_value (27), "Number of grid points, should be odd")
      ("kz", po::value<int > (&K.z)->default_value (27), "Number of grid points, should be odd")
      ("orientation",po::value<std::string > (&ofile)->default_value("conf.gro"), "sample of water orientation")
      ("grid,k",po::value<int > (&kValue), "Number of grid points, should be odd, this will overwrite kx, ky and kz")
      ("rhoMol,r", po::value<double > (&rhoMol)->default_value(32.784), "liquid density")
      ("gr-up", po::value<double > (&rhoMol)->default_value(1.0), "upper cut of g(r)")
      ("rdf-OO",  po::value<std::string > (&fileOO)->default_value ("rdf.oo.1e-2.xvg"), "rdf file of OO")
      ("rdf-OH",  po::value<std::string > (&fileOH)->default_value ("rdf.oh.1e-2.xvg"), "rdf file of OH")
      ("rdf-HH",  po::value<std::string > (&fileHH)->default_value ("rdf.hh.1e-2.xvg"), "rdf file of HH")
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
    printf ("## charge table: %s\n", qfile.c_str());
  }
  printf ("#######################################################\n");

  std::vector<int >  resdindex;
  std::vector<std::string >   resdname;
  std::vector<std::string >   atomname;
  std::vector<int >  atomindex;
  std::vector<std::vector<double > >  posi;
  std::vector<std::vector<double > >  velo;
  std::vector<double >  boxsize;
  GroFileManager::read (ofile, resdindex, resdname, atomname, atomindex, posi, velo, boxsize);
  std::vector<double > charges(posi.size());
  {
    FILE * fptable = fopen(qfile.c_str(), "r");
    if (fptable == NULL){
      std::cerr << "cannot open file " << qfile << std::endl;
      exit (1);
    }
    for (int i = 0; i < int(posi.size()); ++i){
      double tmpvalue;
      int returnvalue = fscanf(fptable, "%lf", &tmpvalue);
      if (returnvalue != 1){
	std::cerr << "wrong format of file " << qfile << std::endl;
	exit(1);
      }
      charges[i] = tmpvalue;
    }
    fclose (fptable);
  }

  DensityProfile_PiecewiseConst dp;

  if (vm.count("charge-table")){
    dp.reinit_xtc (tfile.c_str(), qfile.c_str(), K.x, K.y, K.z, start, end);
  }
  else {
    dp.reinit_xtc (tfile.c_str(), K.x, K.y, K.z, start, end);
  }
  dp.print_avg_x (rfile.c_str());

  ErrorEstimate_SPME_St_H2O_gr eesi;
  eesi.reinit (beta, order, dp, posi, charges, rhoMol, fileOO, fileOH, fileHH, grup);
  eesi.calError (dp, pcharge);
  eesi.print_error (efile.c_str());
  eesi.print_meanf ("meanf.out", dp);

  return 0;
}
