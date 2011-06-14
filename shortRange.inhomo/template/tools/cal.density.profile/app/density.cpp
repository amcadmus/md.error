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
#include "BlockAverage.h"

namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  std::string filename;
  po::options_description desc ("Allow options");  

  float start_t, end_t;
  double h;
  double chl1, chl2, chg1, chg2;
  unsigned nblock;
  
  desc.add_options()
      ("help,h", "print this message")
      ("nblock,n", po::value<unsigned > (&nblock)->default_value (16),  "num of block")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("bin-size,h",  po::value<double > (&h)->default_value (1.),  "bin size")
      ("check-point-liquid-1", po::value<double > (&chl1)->default_value (60.), "liquid check point 1")
      ("check-point-liquid-2", po::value<double > (&chl2)->default_value (90.), "liquid check point 2")
      ("check-point-gas-1", po::value<double > (&chg1)->default_value (40. ), "gas check point 1")
      ("check-point-gas-2", po::value<double > (&chg2)->default_value (110.), "gas check point 2")
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("traj.xtc"), "trajactory file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  DensityProfile_PiecewiseConst dp (filename, h);
  dp.deposit (filename, start_t, end_t);
  
  if (chl1 >= dp.getBox()[0] ||
      chl2 >= dp.getBox()[0] ||
      chg1 >= dp.getBox()[0] ||
      chg2 >= dp.getBox()[0]){
    std::cerr << "check points should be smaller than Lx" << std::endl;
    return 1;
  }

  std::cout << "# start at: " << start_t << std::endl;
  std::cout << "# end   at: " << end_t   << std::endl;
  std::cout << "# liquid from " << chl1 << " to " << chl2 << std::endl;
  std::cout << "# gas    from " << 0    << " to " << chg1 << std::endl;
  std::cout << "# gas    from " << chg2 << " to " << dp.getBox()[0] << std::endl;  
  std::cout << "# profile bin size: " << dp.getBox()[0] / dp.getNx() << std::endl;
  
  BlockAverage avgLiquid, avgGas;
  std::vector<double > denLiquid, denGas;
  denLiquid.reserve (sizeof(double) * dp.getProfile().size());
  double hh = dp.getBox()[0] / dp.getNx();
  for (unsigned ix = 0; ix < dp.getNx(); ++ix){
    double xx = hh * (ix + 0.5);
    if (xx <= chg1) {
      for (unsigned iy = 0; iy < dp.getNy(); ++iy){
	for (unsigned iz = 0; iz < dp.getNz(); ++iz){
	  denGas.push_back (dp.getProfile(ix, iy, iz));
	}
      }
    }
    else if (xx >= chl1 && xx <= chl2){
      for (unsigned iy = 0; iy < dp.getNy(); ++iy){
	for (unsigned iz = 0; iz < dp.getNz(); ++iz){
	  denLiquid.push_back (dp.getProfile(ix, iy, iz));
	}
      }
    }
    else if (xx >= chg2){
      for (unsigned iy = 0; iy < dp.getNy(); ++iy){
	for (unsigned iz = 0; iz < dp.getNz(); ++iz){
	  denGas.push_back (dp.getProfile(ix, iy, iz));
	}
      }
    }
  }

  avgLiquid.processData (denLiquid, nblock);
  avgGas.processData (denGas, nblock);

  printf ("# gas    rho: %.6f ( %.6f )\n", avgGas.getAvg(),    avgGas.getAvgError());
  printf ("# liquid rho: %.6f ( %.6f )\n", avgLiquid.getAvg(), avgLiquid.getAvgError());
  
  return 0;
}

  
