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
  int prec_digit;
  
  desc.add_options()
      ("help,h", "print this message")
      ("nblock,n", po::value<unsigned > (&nblock)->default_value (16),  "num of block")
      ("start,s", po::value<float > (&start_t)->default_value (0.f),  "start time")
      ("end,e",   po::value<float > (&end_t)  ->default_value (0.f),  "end time, 0 is infinity")
      ("bin-size",  po::value<double > (&h)->default_value (1.),  "bin size")
      ("precision,p",  po::value<int > (&prec_digit)->default_value (4),  "number of precision digits")
      ("check-point-liquid-1", po::value<double > (&chl1)->default_value (60.), "liquid check point 1")
      ("check-point-liquid-2", po::value<double > (&chl2)->default_value (90.), "liquid check point 2")
      ("check-point-gas-1", po::value<double > (&chg1)->default_value (40. ), "gas check point 1")
      ("check-point-gas-2", po::value<double > (&chg2)->default_value (110.), "gas check point 2")
      ("file-name,f", po::value<std::string > (&filename)->default_value (std::string("traj.xtc"), "trajactory file name"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help") || vm.count("h")){
    std::cout << desc<< "\n";
    return 0;
  }

  DensityProfile_PiecewiseConst dp (filename, h);
  // dp.deposit (filename, start_t, end_t);
  
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
  double volumeg = dp.getBox()[1] * dp.getBox()[2] * (chg1 + dp.getBox()[0] - chg2);
  double volumel = dp.getBox()[1] * dp.getBox()[2] * (chl2 - chl1);

  XDRFILE * fp = xdrfile_open (filename.c_str(), "r");  
  int step;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * dp.getNatoms());
  float time_prec = 0.001;
  while (read_xtc (fp, dp.getNatoms(), &step, &time, gbox, xx, &prec) == 0){
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
    // calculate rho
    double rhog, rhol;
    rhog = rhol = 0.;
    for (unsigned i = 0; i < unsigned(dp.getNatoms()); ++i) {
      double tmp;
      tmp = xx[i][0];
      if      (xx[i][0] >= dp.getBox()[0]) tmp -= dp.getBox()[0];
      else if (xx[i][0] <  0)              tmp += dp.getBox()[0];
      if (tmp <= chg1 || tmp >= chg2){
	rhog += 1.;
      }
      else if (tmp >= chl1 && tmp <= chl2){
	rhol += 1.;
      }
    }
    if (rhog != 0.) rhog /= volumeg;
    if (rhol != 0.) rhol /= volumel;
    denGas.push_back (rhog);
    denLiquid.push_back (rhol);
  }
  printf ("\n");
  
  avgLiquid.processData (denLiquid, nblock);
  avgGas.processData (denGas, nblock);

  char string[1024];
  sprintf (string, "# gas    rho: %%.%df ( %%.%df ) \n", prec_digit, prec_digit);
  printf  (string, avgGas.getAvg(),    avgGas.getAvgError() * 2.);
  sprintf (string, "# liquid rho: %%.%df ( %%.%df ) \n", prec_digit, prec_digit);
  printf  (string, avgLiquid.getAvg(), avgLiquid.getAvgError() * 2.);
  
  free(xx);
  xdrfile_close(fp);

  return 0;
}

  
