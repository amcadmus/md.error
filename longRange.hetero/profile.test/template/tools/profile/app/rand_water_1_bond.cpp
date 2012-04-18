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
#include "RandomGenerator.h"
#include "GroFileManager.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  double Lx, Ly, Lz;
  double rhoh, rhol;
  double peak;
  double layer;
  std::string cfile, tfile, qfile, wfile;
  int nframe;
  unsigned long seed;
  double charge;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("seed,s", po::value<unsigned long > (&seed)->default_value(1), "random seed")
      ("charge,q", po::value<double > (&charge)->default_value(1.), "abs charge value")
      ("Lx,x", po::value<double > (&Lx)->default_value(30.), "box dim of x")
      ("Ly,y", po::value<double > (&Ly)->default_value(30.), "box dim of y")
      ("Lz,z", po::value<double > (&Lz)->default_value(30.), "box dim of z")
      ("peak,p", po::value<double > (&peak)->default_value(10.), "peak size")
      ("layer,t", po::value<double > (&layer)->default_value(4.), "layer thickness")
      ("rhoh,u", po::value<double > (&rhoh)->default_value(0.50), "higher density")
      ("rhol,l", po::value<double > (&rhol)->default_value(0.05), "lower  density")
      ("nframe,n", po::value<int > (&nframe)->default_value(100), "number of frame in trajctory")
      ("water-conf",  po::value<std::string > (&wfile)->default_value ("water.gro"), "the water conf file")
      ("output-conf",  po::value<std::string > (&cfile)->default_value ("conf.gro"), "the output conf file")
      ("output-charge",po::value<std::string > (&qfile)->default_value ("charge.tab"), "the charge table")
      ("output-traj",  po::value<std::string > (&tfile)->default_value ("traj.xtc"), "the output traj file");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  RandomGenerator_MT19937::init_genrand (seed);

  std::vector<double > yp (6, .5);
  std::vector<double > yn (6, .5);
  std::vector<double > x (6, 20);
  yp[0] = yp[1] = yp[4] = yp[5] = rhol;
  yp[2] = yp[3] = rhoh;
  yn[0] = yn[1] = yn[4] = yn[5] = rhol * 2.;
  yn[2] = yn[3] = rhoh * 2.;
  x[0] = 0.;
  x[1] = 0.5 * (Lx - peak - layer);
  x[2] = 0.5 * (Lx - peak + layer);
  x[3] = x[1] + peak;
  x[4] = x[2] + peak;
  x[5] = Lx;

  SystemDensityFunction sdf;
  sdf.reinit (x, x, yp, yn, Ly, Lz);

  printf ("#######################################################\n");
  printf ("## Lx Ly Lz: %f %f %f\n", Lx, Ly, Lz);
  printf ("## peak:  %f\n", peak);
  printf ("## layer: %f\n", layer);
  printf ("## rhoh: %f\n", rhoh);
  printf ("## rhol: %f\n", rhol);
  printf ("## charge: %f\n", charge);
  printf ("## nframe: %d\n", nframe);
  printf ("## water file: %s\n", wfile.c_str());
  printf ("## n +: %d\n", sdf.get_natom_posi());
  printf ("## n -: %d\n", sdf.get_natom_nega());
  printf ("## n  : %d\n", sdf.get_natom());
  printf ("#######################################################\n");

  std::vector<int > resdindex;
  std::vector<std::string >  resdname;
  std::vector<std::string >  atomname;
  std::vector<int > atomindex;
  std::vector<std::vector<double > > posi;
  std::vector<std::vector<double > > velo;
  std::vector<double > boxsize;

  GroFileManager::read (wfile, resdindex, resdname, atomname, atomindex,
			posi, velo, boxsize);
  
  sdf.genConf_rand_water_gro (cfile.c_str(), posi);
  sdf.genXtc_rand_water (tfile.c_str(), nframe, posi);

  FILE* fp_table = fopen (qfile.c_str(), "w");
  if (fp_table == NULL){
    std::cerr << "cannot open file " << qfile << std::endl;
    return 1;
  }
  for (unsigned i = 0; i < sdf.get_natom_posi(); ++i){
    fprintf (fp_table, "%f\n", charge);
    fprintf (fp_table, "%f\n", - charge * 0.5);
    fprintf (fp_table, "%f\n", - charge * 0.5);
  }
  
  fclose (fp_table);
  
  return 0;
}

