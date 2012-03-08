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

#include "GroFileManager.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  po::options_description desc ("Allow options");
  std::string config, force0, force1, out, meanf, qfile;
  double refh;
  
  desc.add_options()
      ("help,h", "print this message")
      ("config,c", po::value<std::string > (&config)->default_value("conf.gro"), "config file")
      ("charge-table,q", po::value<std::string > (&qfile), "the charge table")
      ("refh", po::value<double > (&refh)->default_value(1.), "bin size")
      ("force-0", po::value<std::string > (&force0)->default_value("force0.out"), "force file 0")
      ("force-1", po::value<std::string > (&force1)->default_value("force1.out"), "force file 1")
      ("error-file,o", po::value<std::string > (&out)->default_value("error.out"), "error file")
      ("mean-force-file,m", po::value<std::string > (&meanf)->default_value("meanf.out"), "error file")
      ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  FILE * fp0, * fp1;
  fp0 = fopen (force0.c_str(), "r");
  fp1 = fopen (force1.c_str(), "r");
  if (fp0 == NULL){
    std::cerr << "cannot open file " << force0 << std::endl;
    return 1;
  }
  if (fp1 == NULL){
    std::cerr << "cannot open file " << force1 << std::endl;
    return 1;
  }
  
  std::vector<int > resdindex, atomindex;
  std::vector<std::string > resdname, atomname;
  std::vector<std::vector<double > > posi, velo;
  std::vector<double > tmpbox;
  GroFileManager::read (config,
			resdindex, resdname, atomname, atomindex,
			posi, velo,
			tmpbox);
  
  int nbin = tmpbox[0] / refh + 0.5;
  double hx = tmpbox[0] / double(nbin);
  // std::vector<int > count(nbin, 0);
  std::vector<int > countp(nbin, 0);

  std::vector<double > dir_value  (nbin, 0.);
  std::vector<double > dir_meanfx (nbin, 0.);
  std::vector<double > dir_meanfy (nbin, 0.);
  std::vector<double > dir_meanfz (nbin, 0.);

  std::vector<double > rec_value  (nbin, 0.);
  std::vector<double > rec_meanfx (nbin, 0.);
  std::vector<double > rec_meanfy (nbin, 0.);
  std::vector<double > rec_meanfz (nbin, 0.);

  std::vector<double > t_value  (nbin, 0.);
  std::vector<double > t_meanfx (nbin, 0.);
  std::vector<double > t_meanfy (nbin, 0.);
  std::vector<double > t_meanfz (nbin, 0.);

  std::vector<double > charge (posi.size(), 1.);
  for (unsigned i = posi.size() / 2; i < posi.size(); ++i){
    charge[i] = -1.;
  }
  if (vm.count("charge-table")){
    FILE * fptable = fopen(qfile.c_str(), "r");
    if (fptable == NULL){
      std::cerr << "cannot open file " << qfile << std::endl;
      exit (1);
    }
    for (unsigned i = 0; i < posi.size(); ++i){
      double tmpvalue;
      int returnvalue = fscanf(fptable, "%lf", &tmpvalue);
      if (returnvalue != 1){
	std::cerr << "wrong format of file " << qfile << std::endl;
	exit(1);
      }
      charge[i] = tmpvalue;
    }
    fclose (fptable);
  }
  
  for (unsigned i = 0; i < posi.size(); ++i){
    if (posi[i][0] >= tmpbox[0]) posi[i][0] -= tmpbox[0];
    else if (posi[i][0] < 0) posi[i][0] += tmpbox[0];
    int idx = posi[i][0] / hx;
    double dir_f0x, dir_f0y, dir_f0z;
    double dir_f1x, dir_f1y, dir_f1z;
    double rec_f0x, rec_f0y, rec_f0z;
    double rec_f1x, rec_f1y, rec_f1z;
    double t_f0x, t_f0y, t_f0z;
    double t_f1x, t_f1y, t_f1z;
    fscanf (fp0, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	    &dir_f0x, &dir_f0y, &dir_f0z,
	    &rec_f0x, &rec_f0y, &rec_f0z,
	    &t_f0x, &t_f0y, &t_f0z
	);
    fscanf (fp1, "%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	    &dir_f1x, &dir_f1y, &dir_f1z,
	    &rec_f1x, &rec_f1y, &rec_f1z,
	    &t_f1x, &t_f1y, &t_f1z
	);
    dir_f0x -= dir_f1x;
    dir_f0y -= dir_f1y;
    dir_f0z -= dir_f1z;
    rec_f0x -= rec_f1x;
    rec_f0y -= rec_f1y;
    rec_f0z -= rec_f1z;
    t_f0x -= t_f1x;
    t_f0y -= t_f1y;
    t_f0z -= t_f1z;
    if (charge[i] > 0){
      dir_meanfx[idx] += dir_f0x;
      dir_meanfy[idx] += dir_f0y;
      dir_meanfz[idx] += dir_f0z;
      rec_meanfx[idx] += rec_f0x;
      rec_meanfy[idx] += rec_f0y;
      rec_meanfz[idx] += rec_f0z;
      t_meanfx[idx] += t_f0x;
      t_meanfy[idx] += t_f0y;
      t_meanfz[idx] += t_f0z;
      dir_value[idx] += dir_f0x*dir_f0x + dir_f0y*dir_f0y + dir_f0z*dir_f0z;
      rec_value[idx] += rec_f0x*rec_f0x + rec_f0y*rec_f0y + rec_f0z*rec_f0z;
      t_value[idx] += t_f0x*t_f0x + t_f0y*t_f0y + t_f0z*t_f0z;
      countp[idx] ++;
    }
    // else {
    //   meanfx[idx] -= f0x;
    //   meanfy[idx] -= f0y;
    //   meanfz[idx] -= f0z;
    //   countp[idx] ++;
    // }
    // count[idx] += 1;
  }
  for (int i = 0; i < nbin; ++i){
    if (countp[i] != 0){
      dir_meanfx[i] /= countp[i];
      dir_meanfy[i] /= countp[i];
      dir_meanfz[i] /= countp[i];
      rec_meanfx[i] /= countp[i];
      rec_meanfy[i] /= countp[i];
      rec_meanfz[i] /= countp[i];
      t_meanfx[i] /= countp[i];
      t_meanfy[i] /= countp[i];
      t_meanfz[i] /= countp[i];
      dir_value[i] = sqrt (dir_value[i] / countp[i]);
      rec_value[i] = sqrt (rec_value[i] / countp[i]);
      t_value[i] = sqrt (t_value[i] / countp[i]);
    }
    // if (count[i] != 0){
    // }
  }
  
  fclose (fp0);
  fclose (fp1);
  
  FILE * fp = fopen(out.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << out << std::endl;
    return 1;
  }
  fprintf(fp, "# 1  2      3      4      5      \n");
  fprintf(fp, "# x  dirEF  recEF  sumEF  #Sample\n");	  
  for (int i = 0; i < nbin; ++i){
    fprintf(fp, "%f %e %e %e %d\n",
	    (i+0.5) * refh,
	    dir_value[i],
	    rec_value[i],
	    t_value[i],
	    countp[i]
	);
  }
  fclose (fp);
  fp = fopen(meanf.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << out << std::endl;
    return 1;
  }
  fprintf(fp, "# 1  2-4    5-7    8-10   11     \n");
  fprintf(fp, "# x  dirEF  recEF  sumEF  #Sample\n");	  
  for (int i = 0; i < nbin; ++i){
    fprintf(fp, "%f \t %e %e %e \t %e %e %e \t %e %e %e \t %d\n",
	    (i+0.5) * refh,
	    dir_meanfx[i],
	    dir_meanfy[i],
	    dir_meanfz[i],
	    rec_meanfx[i],
	    rec_meanfy[i],
	    rec_meanfz[i],
	    t_meanfx[i],
	    t_meanfy[i],
	    t_meanfz[i],
	    countp[i]
	);
  }
  fclose (fp);
}


