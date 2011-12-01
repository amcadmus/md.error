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
#include "Interpolation.h"
#include "StringSplit.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void readMFile (const std::string & file,
		const int & nx,
		const int & vx,
		const int & vy,
		const int & vz,
		std::vector<double > & xx,
		std::vector<double > & vvx,
		std::vector<double > & vvy,
		std::vector<double > & vvz)
{
  FILE * fp = fopen (file.c_str(), "r");
  if (NULL == fp){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }
  char line[2048];
  xx.clear();
  vvx.clear();
  vvy.clear();
  vvz.clear();
  while (fgets (line, 2048, fp) != NULL){
    if (line[0] == '#') continue;
    std::vector<std::string > words;
    StringOperation::split (line, words);
    if ((nx-1) >= int(words.size()) || (vx-1) >= int(words.size()) ||
	(vy-1) >= int(words.size()) || (vz-1) >= int(words.size()) ){
      std::cerr << "file fomat of " << file << " is wrong" << std::endl;
      exit(1);
    }
    xx. push_back(atof(words[(nx-1)].c_str()));
    vvx.push_back(atof(words[(vx-1)].c_str()));
    vvy.push_back(atof(words[(vy-1)].c_str()));
    vvz.push_back(atof(words[(vz-1)].c_str()));
  }
  fclose (fp);
}

	
void readEFile (const std::string & file,
		const int & nx,
		const int & vx,
		std::vector<double > & xx,
		std::vector<double > & vv)
{
  FILE * fp = fopen (file.c_str(), "r");
  if (NULL == fp){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }
  char line[2048];
  xx.clear();
  vv.clear();
  while (fgets (line, 2048, fp) != NULL){
    if (line[0] == '#') continue;
    std::vector<std::string > words;
    StringOperation::split (line, words);
    if ((nx-1) >= int(words.size()) || (vx-1) >= int(words.size()) ){
      std::cerr << "file fomat of " << file << " is wrong" << std::endl;
      exit(1);
    }
    xx.push_back(atof(words[(nx-1)].c_str()));
    vv.push_back(atof(words[(vx-1)].c_str()));
  }
  fclose (fp);
}




int main(int argc, char * argv[])
{
  std::string ofile, ffile;
  std::string dir_meanf, rec_meanf, dir_error, rec_error;
  double refh;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("refh,n", po::value<double > (&refh)->default_value(1.), "bin size")
      ("dir-meanf", po::value<std::string > (&dir_meanf)->default_value("dir.meanf.out"), "mean force of the dir part")
      ("rec-meanf", po::value<std::string > (&rec_meanf)->default_value("rec.meanf.out"), "mean force of the rec part")
      ("dir-error", po::value<std::string > (&dir_error)->default_value("dir.error.out"), "RMS error of the dir part")
      ("rec-error", po::value<std::string > (&rec_error)->default_value("rec.error.out"), "RMS error of the rec part")
      ("output-error",  po::value<std::string > (&ofile)->default_value ("error.out"), "the output error of the system")
      ("output-meanf",  po::value<std::string > (&ffile)->default_value ("meanf.out"), "the output error of the system")
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::vector<double > dmx, rmx;
  std::vector<double > dmvx, dmvy, dmvz;
  std::vector<double > rmvx, rmvy, rmvz;
  std::vector<double > de, re;
  
  readMFile (dir_meanf, 1, 3, 4, 5, dmx, dmvx, dmvy, dmvz);
  readMFile (rec_meanf, 1, 3, 4, 5, rmx, rmvx, rmvy, rmvz);
  readEFile (dir_error, 1, 2, dmx, de);
  readEFile (rec_error, 1, 2, rmx, re);

  PiecewisePoly dmeanfx, dmeanfy, dmeanfz;
  PiecewisePoly rmeanfx, rmeanfy, rmeanfz;
  PiecewisePoly derror,  rerror;
  
  Interpolation::piecewiseLinearPeriodic (dmx, dmvx, dmeanfx);
  Interpolation::piecewiseLinearPeriodic (dmx, dmvy, dmeanfy);
  Interpolation::piecewiseLinearPeriodic (dmx, dmvz, dmeanfz);
  Interpolation::piecewiseLinearPeriodic (rmx, rmvx, rmeanfx);
  Interpolation::piecewiseLinearPeriodic (rmx, rmvy, rmeanfy);
  Interpolation::piecewiseLinearPeriodic (rmx, rmvz, rmeanfz);

  Interpolation::piecewiseLinearPeriodic (dmx, de, derror);
  Interpolation::piecewiseLinearPeriodic (rmx, re, rerror);

  // FILE * fp = fopen("test.out", "w");
  // double x0 = 0;
  // double x1=  30.;
  // double nn = 1000;
  // double h = (x1 - x0) / double(nn);
  // for (int i = 0; i <= nn; ++i){
  //   fprintf (fp, "%f %f %f\n",
  // 	     x0 + i * h,
  // 	     dmeanfx.value(x0 + i * h),
  // 	     rmeanfx.value(x0 + i * h));
  // }
  // fclose (fp);

  double box = dmx.back() + ( 0.5 * (dmx[1] - dmx[0]) );
  int nn = unsigned (box / refh);
  double hh = box / double (nn);
  
  FILE * fp = fopen (ofile.c_str(), "w");
  fprintf (fp, "# x  Error  E_dir  E_rec  2*mF_dir*mF_rec\n");
  for (int i = 0; i < nn; ++i){
    double xx = (i+0.5) * hh;
    double tmp1 = derror.value (xx);
    tmp1 = tmp1 * tmp1;
    double tmp2 = rerror.value (xx);
    tmp2 = tmp2 * tmp2;
    double tmp3 = 2. * (
	dmeanfx.value(xx) * rmeanfx.value(xx) + 
	dmeanfy.value(xx) * rmeanfy.value(xx) + 
	dmeanfz.value(xx) * rmeanfz.value(xx) );
	
    fprintf (fp, "%f %e    %e %e %e\n",
	     xx, sqrt (tmp1 + tmp2 + tmp3),
	     sqrt(tmp1), sqrt(tmp2), tmp3);
  }
  fclose (fp);

  fp = fopen (ffile.c_str(), "w");
  fprintf (fp, "# x  mFx  mFy  mFz,   \n");
  fprintf(fp, "# 1  2-4    5-7    8-10   11     \n");
  fprintf(fp, "# x  sumEF  dirEF  recEF  #Sample\n");	  
  for (int i = 0; i < nn; ++i){
    double xx = (i+0.5) * hh;
    double fx = dmeanfx.value(xx) + rmeanfx.value(xx);
    double fy = dmeanfy.value(xx) + rmeanfy.value(xx);
    double fz = dmeanfz.value(xx) + rmeanfz.value(xx);
    fprintf (fp, "%f     %e %e %e     %e %e %e     %e %e %e\n",
	     xx, 
	     fx, fy, fz,
	     dmeanfx.value(xx), dmeanfy.value(xx), dmeanfz.value(xx),
	     rmeanfx.value(xx), rmeanfy.value(xx), rmeanfz.value(xx)
	);    
  }
  fclose (fp);
  
  return 0;
}
