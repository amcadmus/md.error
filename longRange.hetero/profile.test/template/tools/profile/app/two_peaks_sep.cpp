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
  double layer;
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
      ("layer,t", po::value<double > (&layer)->default_value(4.), "layer thickness")
      ("rhoh", po::value<double > (&rhoh)->default_value(0.50), "charge density")
      ("rhol", po::value<double > (&rhol)->default_value(0.05), "charge density")
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
  printf ("## print two peaks seperating the box\n");
  printf ("## Lx Ly Lz: %f %f %f\n", Lx, Ly, Lz);
  printf ("## layer: %f\n", layer);
  printf ("## rhoh: %f\n", rhoh);
  printf ("## rhol: %f\n", rhol);
  printf ("## nframe: %d\n", nframe);
  printf ("#######################################################\n");

  std::vector<double > yp (6, .5);
  std::vector<double > yn (6, .0);
  std::vector<double > xp (6, 20);
  std::vector<double > xn (6, 20);

  xp[0] = 0.0;
  xp[1] = 0.25 * Lx - 0.50 * layer;
  xp[2] = 0.25 * Lx + 0.50 * layer;
  xp[3] = 0.75 * Lx - 0.50 * layer;
  xp[4] = 0.75 * Lx + 0.50 * layer;
  xp[5] = Lx;

  xn = xp;

  yp[0] = rhol;
  yp[1] = rhol;
  yp[2] = rhoh;
  yp[3] = rhoh;
  yp[4] = rhol;
  yp[5] = rhol;

  yn[0] = rhoh;
  yn[1] = rhoh;
  yn[2] = rhol;
  yn[3] = rhol;
  yn[4] = rhoh;
  yn[5] = rhoh;
  
  SystemDensityFunction sdf;
  sdf.reinit (xp, xn, yp, yn, Ly, Lz);
  sdf.genConf_gro (cfile.c_str());
  sdf.genXtc (tfile.c_str(), nframe);
  return 0;
}

