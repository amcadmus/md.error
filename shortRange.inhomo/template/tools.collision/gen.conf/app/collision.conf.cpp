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
#include "GroFileManager.h"

namespace po = boost::program_options;

void analyze (std::vector<std::vector<double > > & posi,
	      std::vector<std::vector<double > > & velo,
	      double & rmax,
	      const double & rthreadhold,
	      std::vector<double > & center)
{
  restart_line:
  
  center.resize (3);
  center[0] = center[1] = center[2] = 0.;
  for (unsigned i = 0; i < posi.size(); ++i){
    center[0] += posi[i][0];
    center[1] += posi[i][1];
    center[2] += posi[i][2];
  }
  center[0] /= double (posi.size());
  center[1] /= double (posi.size());
  center[2] /= double (posi.size());

  std::vector<std::vector<double > >::iterator itp = posi.begin();
  std::vector<std::vector<double > >::iterator itv = velo.begin();
  rmax = 0;
  for (; itp != posi.end(); ++itp, ++itv){
    double dx = itp->operator[](0) - center[0];
    double dy = itp->operator[](1) - center[1];
    double dz = itp->operator[](2) - center[2];
    double diff = sqrt (dx*dx + dy*dy + dz*dz);
    if (diff > rthreadhold){
      std::cout << "# kick atom, the diff is " << diff << std::endl;
      posi.erase (itp);
      velo.erase (itv);
      goto restart_line;
    }
    if (diff > rmax) rmax = diff;
  }
}

  
int main (int argc, char * argv[])
{
  double Lx, Ly, Lz;
  double xx, uu;
  
  std::string inputfile;
  std::string filename;
  po::options_description desc ("Allow options");  

  desc.add_options()
      ("help,h", "print this message")
      ("Lx", po::value<double > (&Lx)->default_value (150),  "box size x.")
      ("Ly", po::value<double > (&Ly)->default_value (27*2),  "box size y.")
      ("Lz", po::value<double > (&Lz)->default_value (27*3),  "box size z.")
      ("x-value,x", po::value<double > (&xx)->default_value (0.),  "value of x.")
      ("u-value,u", po::value<double > (&uu)->default_value (1.58),"value of u.")
      ("ball-conf",	po::value<std::string > (&inputfile)->	default_value (std::string("ball.init.gro"), "inital ball conf."))
      ("file-name,f",	po::value<std::string > (&filename)->	default_value (std::string("confout.gro"), "file to output."));
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::vector<std::vector<double > > posi, velo;
  std::vector<double > zero (3, 0.);
  std::vector<std::string > resdname, atomname;
  std::vector<int > resdindex, atomindex;
  std::vector<double > box(3);

  GroFileManager::read (inputfile,
			resdindex, resdname,
			atomname, atomindex,
			posi, velo, box);
  std::vector<double > center;
  double rmax;
  double rthreadhold = 30;
  analyze (posi, velo, rmax, rthreadhold, center);
  double dmax = 2.*rmax;
  std::cout << "# rmax is " << rmax << std::endl;
  std::cout << "# dmax is " << dmax << std::endl;
  std::cout << "# x value is " << xx << std::endl;
  std::cout << "# u value is " << uu << std::endl;
  
  double tmpdmax = 28;
  Ly = 2. * tmpdmax;
  Lz = 3. * tmpdmax;
  double XX = xx * tmpdmax;
  std::vector<double > newcenter0(3), newcenter1(3);
  newcenter0[0] = tmpdmax;
  newcenter0[1] = 0.5 * Ly;
  newcenter0[2] = 0.5 * Lz + 0.5 * XX;
  newcenter1[0] = Lx - tmpdmax;
  newcenter1[1] = 0.5 * Ly;
  newcenter1[2] = 0.5 * Lz - 0.5 * XX;
  
  std::vector<std::vector<double > > posi_1, velo_1;
  std::vector<double > shift (3);
  shift[0] = newcenter0[0] - center[0];
  shift[1] = newcenter0[1] - center[1];
  shift[2] = newcenter0[2] - center[2];
  for (unsigned i = 0; i < posi.size(); ++i) {
    std::vector<double > tmpp (posi[i]);
    tmpp[0] += shift[0];
    tmpp[1] += shift[1];
    tmpp[2] += shift[2];
    posi_1.push_back (tmpp);
    std::vector<double > tmpv (velo[i]);
    tmpv[0] += 0.5 * uu;
    velo_1.push_back (tmpv);
  }
  shift[0] = newcenter1[0] - center[0];
  shift[1] = newcenter1[1] - center[1];
  shift[2] = newcenter1[2] - center[2];
  for (unsigned i = 0; i < posi.size(); ++i) {
    std::vector<double > tmpp (posi[i]);
    tmpp[0] += shift[0];
    tmpp[1] += shift[1];
    tmpp[2] += shift[2];
    posi_1.push_back (tmpp);
    std::vector<double > tmpv (velo[i]);
    tmpv[0] -= 0.5 * uu;
    velo_1.push_back (tmpv);
  }
  
  resdname.resize (posi_1.size());
  resdindex.resize (posi_1.size());

  unsigned count = 1;
  for (unsigned i = 0; i < posi_1.size()/2; ++i){
    resdname[i] = (std::string("lja"));
    resdindex[i] = count;
    count ++;
  }
  for (unsigned i = posi_1.size()/2; i < posi_1.size(); ++i){
    resdname[i] = (std::string("ljb"));
    resdindex[i] = count;
    count ++;
  }

  box[0] = Lx;
  box[1] = Ly;
  box[2] = Lz;
  
  GroFileManager::write (filename,
			 resdindex, resdname,
			 resdname, resdindex,
			 posi_1, velo_1, box);
  return 0;
}

