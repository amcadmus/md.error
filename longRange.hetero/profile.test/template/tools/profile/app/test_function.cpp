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


int main(int argc, char * argv[])
{
  std::vector<double > y (6, .5);
  std::vector<double > x (6, 20);
  y[0] = 0.1;
  y[1] = 0.1;
  y[2] = 0.5;
  y[3] = 0.5;
  y[4] = 0.1;
  y[5] = 0.1;
  x[0] = 0.;
  x[1] = 8.;
  x[2] = 12.;
  x[3] = 18.;
  x[4] = 22.;
  x[5] = 30.;

  SystemDensityFunction sdf;
  sdf.reinit (x, x, y, y, 30., 30.);
  sdf.genConf_gro ("conf.gro");
  sdf.genXtc ("traj.xtc", 100);

  // DensityFunction df;
  // df.reinit (x, y);

  // int n = 1000;
  // double refh = df.Lx() / double(n);
  // FILE * fp = fopen ("f.out", "w");
  // for (int i = 0; i < n; ++i){
  //   double x = 0 + refh * i;
  //   fprintf (fp, "%f %f\n", x, df(x));
  // }
  
  // fclose (fp);
}

