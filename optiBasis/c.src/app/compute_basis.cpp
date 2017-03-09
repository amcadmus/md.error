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

// libs
#include "TableFileLoader.h"
#include "HatSymmHermiteBase.h"
#include "IkError.h"
#include "LossFunc.h"
#include "CardinalBSpline.h"
#include "SimulationRegion.h"

// third party
#include <boost/program_options.hpp>
#include <dlib/optimization.h>

namespace po = boost::program_options;

LossFunc<IkError<HatSymmHermiteBase> > * global_p_loss;


double func_value (const column_vector& m) 
{
  return global_p_loss->value (m);
}

const column_vector func_deriv (const column_vector& m) 
{
  return global_p_loss->deriv (m);
}

void print_result (const char * file,
		   const int & CC,
		   const column_vector & xx,
		   const string & head)
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    cerr << "cannot open file " << file << endl;
    exit(1);
  }
  fprintf (fp, "# %s\n", head.c_str());
  long size = xx.size()/2;
  long nbins = size + 1;
  double hh = CC / double(nbins);
  fprintf (fp, "%f %.16e %.16e\n", 0., 1., 0.);
  for (long ii = 0; ii < size; ++ii){
    fprintf (fp, "%f %.16e %.16e\n", (ii+1) * hh, xx(ii), xx(ii+size));    
  }
  fprintf (fp, "%f %.16e %.16e\n", double(CC), 0., 0.);
  fclose (fp);
}

int main(int argc, char * argv[])
{
  string ifile, ofile;
  string b_style;
  int nbins, CC, sKK, l_cut;
  double beta, LL, tol;
  bool is_vbs = false;
  int numb_threads;
  
  po::options_description desc ("# Allow options");
  desc.add_options()
      ("help,h", "Print this message and exit.")
      ("verbose,v", "Being loud and noisy.")
      ("input,f",	po::value<string > (&ifile), "The input guess. If not set use b-spline as initial guess.")
      ("b-style,B",	po::value<string > (&b_style)->default_value ("SPME"), "The style of B-function, can be SPME or PPPM.")
      ("numb-bins,n",	po::value<int > (&nbins)->default_value (10), "The number of dicretizing bins for basis.")
      ("cut-off,c",	po::value<int > (&CC)->default_value (2), "The cut-off of basis.")
      ("beta,b",	po::value<double > (&beta)->default_value (3.0), "The splitting parameter.")
      ("numb-grid,k",	po::value<int > (&sKK)->default_value (32), "The number of grid points.")
      ("box-size,l",	po::value<double > (&LL)->default_value (3.72412), "The box size.")
      ("tolerence,t",	po::value<double > (&tol)->default_value (1e-6), "The relative tol for opti.")
      ("l-cut",		po::value<int > (&l_cut), "the cut-off of l sum. Guess if not set")
      ("numb-threads,T",po::value<int > (&numb_threads)->default_value (1), "the number of threads")
      ("output,o",	po::value<string > (&ofile)->default_value ("basis.out"), "The output basis.");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    cout << desc<< "\n";
    return 0;
  }
  if (vm.count("verbose")){
    is_vbs = true;
  }
  if (!vm.count("l-cut")){
    l_cut = int (4* nbins / CC);
    if (is_vbs) cout << "# guessed l_cut " << l_cut << endl;
  }

  int ib_style = -1;
  if (b_style == "SPME") {
    ib_style = B_STYLE::SPME;
  }
  else if (b_style == "PPPM") {
    ib_style = B_STYLE::PPPM;    
  }
  else {
    cerr << "unknow b-style " << b_style << endl;
    return 1;
  }

  vector<double > vi (nbins-1), di(nbins-1);
  double hh = CC / double(nbins);
  int order = CC * 2;
  if (!vm.count("input")) {
    Basis<double> *pcbs[15];
    pcbs[ 2] = new CardinalBSpline<double, 2> ();
    pcbs[ 4] = new CardinalBSpline<double, 4> ();
    pcbs[ 6] = new CardinalBSpline<double, 6> ();
    pcbs[ 8] = new CardinalBSpline<double, 8> ();
    pcbs[10] = new CardinalBSpline<double,10> ();
    pcbs[12] = new CardinalBSpline<double,12> ();
    pcbs[14] = new CardinalBSpline<double,14> ();
    double scale =  1./ pcbs[order]->value (0);
    for (int ii = 0; ii < nbins-1; ++ii){
      double xx = (ii+1) * hh;
      vi[ii] = pcbs[order]->value (xx) * scale;
      di[ii] = pcbs[order]->derivative (xx) * scale;
    }
    if (is_vbs) cout << "# use b-spline basis as initial guess" << endl;
  }
  else {
    TableFileLoader tbl (ifile.c_str());
    std::vector<unsigned> cols (3);
    cols[0] = 1;
    cols[1] = 2;
    cols[2] = 3;
    tbl.setColumns (cols);
    std::vector<std::vector<double > > data;
    tbl.loadAll (data);
    vi = data[1];
    di = data[2];
    vi.erase (vi.begin());
    vi.erase (vi.end()-1);
    di.erase (di.begin());
    di.erase (di.end()-1);
    assert (int(vi.size()) == nbins-1);
    assert (int(di.size()) == nbins-1);
    if (is_vbs) cout << "# use " << ifile << " as initial guess" << endl;
  }

  double box[9] = {0};
  box[0] = box[4] = box[8] = LL;
  SimulationRegion<double > region;
  region.reinitBox (box);  
  double q2 = 33.456 * LL*LL*LL * (-0.8476 * -0.8476 + 0.4238 * 0.4238 * 2);
  int natoms = 33.456 * LL*LL*LL * 3;
  vector<int> KK(3, sKK);
      
  LossFunc<IkError<HatSymmHermiteBase> > lf (CC, nbins, beta, KK, q2, natoms, region, l_cut, ib_style, numb_threads);
  global_p_loss = &lf;
  
  column_vector xx (vi.size() * 2);
  for (unsigned ii = 0; ii < vi.size(); ++ii){
    xx(ii) = vi[ii];
    xx(ii+vi.size()) = di[ii];
  }

  // print_result ("init.out", CC, xx, "");
  double initv = func_value (xx);
  if (is_vbs) cout << "# init value " << initv << endl;

  dlib::objective_delta_stop_strategy stop_strategy (tol * initv);
  if (is_vbs) stop_strategy.be_verbose();

  double vmin = 
      dlib::find_min(dlib::bfgs_search_strategy(),
		     stop_strategy,
		     func_value, func_deriv, xx, -1);

  if (is_vbs) cout << "# final value " << vmin << endl;
  
  char head[1024];
  sprintf (head, "%f %f %e", beta, LL/double(sKK), vmin);
  print_result (ofile.c_str(), CC, xx, head);

  return 0;
}
