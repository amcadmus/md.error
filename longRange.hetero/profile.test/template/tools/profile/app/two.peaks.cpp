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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  double Lx, Ly, Lz;
  double rhoh, rhol;
  double peak;
  double layer;
  double dist;
  std::string cfile, tfile;
  int nframe;
  unsigned long seed;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("seed,s", po::value<unsigned long > (&seed)->default_value(1), "random seed")
      ("Lx,x", po::value<double > (&Lx)->default_value(30.), "box dim of x")
      ("Ly,y", po::value<double > (&Ly)->default_value(30.), "box dim of y")
      ("Lz,z", po::value<double > (&Lz)->default_value(30.), "box dim of z")
      ("peak,p", po::value<double > (&peak)->default_value(10.), "peak size")
      ("layer,t", po::value<double > (&layer)->default_value(4.), "layer thickness")
      ("dist,d", po::value<double > (&dist)->default_value(5.), "dist between two peaks")
      ("rhoh,u", po::value<double > (&rhoh)->default_value(0.50), "higher density")
      ("rhol,l", po::value<double > (&rhol)->default_value(0.05), "lower  density")
      ("nframe,n", po::value<int > (&nframe)->default_value(100), "number of frame in trajctory")
      ("output-conf",  po::value<std::string > (&cfile)->default_value ("conf.gro"), "the output conf file")
      ("output-traj",  po::value<std::string > (&tfile)->default_value ("traj.xtc"), "the output traj file");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  RandomGenerator_MT19937::init_genrand (seed);

  printf ("#######################################################\n");
  printf ("## Lx Ly Lz: %f %f %f\n", Lx, Ly, Lz);
  printf ("## peak:  %f\n", peak);
  printf ("## layer: %f\n", layer);
  printf ("## rhoh: %f\n", rhoh);
  printf ("## rhol: %f\n", rhol);
  printf ("## nframe: %d\n", nframe);
  printf ("#######################################################\n");

  std::vector<double > y (6, .5);
  std::vector<double > yn (6, .0);
  std::vector<double > xp (6, 20);
  std::vector<double > xn (6, 20);
  y[0] = y[1] = y[4] = y[5] = rhol;
  y[2] = y[3] = rhoh;
  xp[0] = 0.;
  xp[1] = 0.5 * (Lx - peak - layer) - 0.5 * dist;
  xp[2] = 0.5 * (Lx - peak + layer) - 0.5 * dist;
  xp[3] = xp[1] + peak;
  xp[4] = xp[2] + peak;
  xp[5] = Lx;
  xn[0] = 0.;
  xn[1] = 0.5 * (Lx - peak - layer) + 0.5 * dist;
  xn[2] = 0.5 * (Lx - peak + layer) + 0.5 * dist;
  xn[3] = xn[1] + peak;
  xn[4] = xn[2] + peak;
  xn[5] = Lx;

  SystemDensityFunction sdf;
  sdf.reinit (xp, xn, y, y, Ly, Lz);
  sdf.genConf_gro (cfile.c_str());
  sdf.genXtc (tfile.c_str(), nframe);
  
  return 0;
}

