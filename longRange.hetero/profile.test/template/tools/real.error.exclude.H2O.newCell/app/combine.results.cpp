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
  std::string record;
  std::string out, meanf;
  double refh=0.;
  
  desc.add_options()
      ("help,h", "print this message")
      ("record-file,r", po::value<std::string > (&record)->default_value("ewald.record"), "file of filenames to be analyzed")
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

  char errorFileName[1024];
  char meanfFileName[1024];

  FILE * fp = fopen (record.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << record << std::endl;
    return 1;
  }

  std::vector<int > count;
  std::vector<int > countp;

  std::vector<double > dir_value  ;
  std::vector<double > dir_meanfx ;
  std::vector<double > dir_meanfy ;
  std::vector<double > dir_meanfz ;

  std::vector<double > rec_value  ;
  std::vector<double > rec_meanfx ;
  std::vector<double > rec_meanfy ;
  std::vector<double > rec_meanfz ;

  std::vector<double > t_value  ;
  std::vector<double > t_meanfx ;
  std::vector<double > t_meanfy ;
  std::vector<double > t_meanfz ;

  int countFile = 0;
  int nbin = 0;
  
  while (fscanf (fp, "%s%s\n", errorFileName, meanfFileName) != EOF){
    FILE * fp0, * fp1;
    char tmpline[1024];

    fp0 = fopen (errorFileName, "r");
    fp1 = fopen (meanfFileName, "r");
    if (fp0 == NULL){
      std::cerr << "cannot open file " << errorFileName << std::endl;
      return 1;
    }
    if (fp1 == NULL){
      std::cerr << "cannot open file " << meanfFileName << std::endl;
      return 1;
    }

    if (0 == countFile){
      countFile ++;
      while (fgets(tmpline, 1024, fp0) != NULL) {
	nbin ++;
      }
      nbin -= 2;
      fclose (fp0);
      fp0 = fopen (errorFileName, "r");
      if (fp0 == NULL){
	std::cerr << "cannot open file " << errorFileName << std::endl;
	return 1;
      }
      count.resize (nbin, 0);
      countp.resize (nbin, 0);
      dir_value.resize (nbin, 0.);
      rec_value.resize (nbin, 0.);
      t_value.resize (nbin, 0.);
      dir_meanfx.resize (nbin, 0.);
      dir_meanfy.resize (nbin, 0.);
      dir_meanfz.resize (nbin, 0.);
      rec_meanfx.resize (nbin, 0.);
      rec_meanfy.resize (nbin, 0.);
      rec_meanfz.resize (nbin, 0.);
      t_meanfx.resize (nbin, 0.);
      t_meanfy.resize (nbin, 0.);
      t_meanfz.resize (nbin, 0.);
    }

    fgets(tmpline, 1024, fp0);
    fgets(tmpline, 1024, fp0);
    fgets(tmpline, 1024, fp1);
    fgets(tmpline, 1024, fp1);

    double mdir_value, mrec_value, mt_value;
    double mdir_meanfx, mdir_meanfy, mdir_meanfz;
    double mrec_meanfx, mrec_meanfy, mrec_meanfz;
    double mt_meanfx, mt_meanfy, mt_meanfz;
    double oldx, tmpx;
    int tmpn;
    oldx = tmpx = 0.;
    
    for (int i = 0; i < nbin; ++i){
      oldx = tmpx;
      fscanf (fp0, "%lf %lf%lf%lf %d",
	      &tmpx,
	      &mdir_value, &mrec_value, &mt_value,
	      &tmpn);
      dir_value[i] += mdir_value * mdir_value * tmpn;
      rec_value[i] += mrec_value * mrec_value * tmpn;
      t_value[i] += mt_value * mt_value * tmpn;
      count[i] += tmpn;

      fscanf (fp1, "%lf %lf%lf%lf %lf%lf%lf %lf%lf%lf %d",
	      &tmpx,
	      &mdir_meanfx, &mdir_meanfy, &mdir_meanfz,
	      &mrec_meanfx, &mrec_meanfy, &mrec_meanfz,
	      &mt_meanfx, &mt_meanfy, &mt_meanfz,
	      &tmpn);
      dir_meanfx[i] += mdir_meanfx * tmpn;
      dir_meanfy[i] += mdir_meanfy * tmpn;
      dir_meanfz[i] += mdir_meanfz * tmpn;
      rec_meanfx[i] += mrec_meanfx * tmpn;
      rec_meanfy[i] += mrec_meanfy * tmpn;
      rec_meanfz[i] += mrec_meanfz * tmpn;
      t_meanfx[i] += mt_meanfx * tmpn;
      t_meanfy[i] += mt_meanfy * tmpn;
      t_meanfz[i] += mt_meanfz * tmpn;
      countp[i] += tmpn;
    }

    refh = tmpx - oldx;
    fclose (fp0);
    fclose (fp1);
  }
  fclose (fp);

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
    }
    if (count[i] != 0){
      dir_value[i] = sqrt (dir_value[i] / count[i]);
      rec_value[i] = sqrt (rec_value[i] / count[i]);
      t_value[i] = sqrt (t_value[i] / count[i]);
    }
  }
  
  
  fp = fopen(out.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << out << std::endl;
    return 1;
  }
  fprintf(fp, "# 1  2      3      4      5      \n");
  fprintf(fp, "# x  sumEF  dirEF  recEF  #Sample\n");	  
  for (int i = 0; i < nbin; ++i){
    fprintf(fp, "%f %e %e %e %d\n",
	    (i+0.5) * refh,
	    t_value[i],
	    dir_value[i],
	    rec_value[i],
	    count[i]
	);
  }
  fclose (fp);
  fp = fopen(meanf.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << out << std::endl;
    return 1;
  }
  fprintf(fp, "# 1  2-4    5-7    8-10   11     \n");
  fprintf(fp, "# x  sumEF  dirEF  recEF  #Sample\n");	  
  for (int i = 0; i < nbin; ++i){
    fprintf(fp, "%f \t %e %e %e \t %e %e %e \t %e %e %e \t %d\n",
	    (i+0.5) * refh,
	    t_meanfx[i],
	    t_meanfy[i],
	    t_meanfz[i],
	    dir_meanfx[i],
	    dir_meanfy[i],
	    dir_meanfz[i],
	    rec_meanfx[i],
	    rec_meanfy[i],
	    rec_meanfz[i],
	    countp[i]
	);
  }
  fclose (fp);
}


