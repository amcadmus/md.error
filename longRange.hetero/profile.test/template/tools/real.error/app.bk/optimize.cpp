#include "Particle.h"
#include "ToolBox.h"
#include "ShortRangePairFinder.h"
#include "BoxGeometry.h"
#include "PairInteraction.h"
#include "VectorOperation.h"
#include "ElectrostaticInteraction.h"
#include "Integral1D.h"
#include "Polynominal.h"
#include "Stopwatch.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

static double unitScale = 138.935485;

int main (int argc, char * argv[]) 
{
//   Stopwatch sw;
//   double a = 1000;
//   double b = 19;
//   double c;
//   sw.start();
//   for (double i = 0; i < 1e8; i = i + 1.) c = a * b;
//   sw.stop();  
//   std::cout << "user time " << sw.user() << std::endl;
//   std::cout << "system time " << sw.system() << std::endl;
//   std::cout << "real time " << sw.real() << std::endl;

// //   sw.start();
//   for (double i = 0; i < 1e8; i = i + 1.) c = a * b;
//   sw.stop();  
//   std::cout << "user time " << sw.user() << std::endl;
//   std::cout << "system time " << sw.system() << std::endl;
//   std::cout << "real time " << sw.real() << std::endl;
  

//   return 0;
  
  unsigned seed;
  double beta;
  unsigned ksize;
  double rcut;
  double boxsize;
  unsigned npart;
  std::string method;
  unsigned order;
  po::options_description desc ("Allow options");
  std::string config;
  std::string forceFile;
  int testInteger;
  double required;
  
  desc.add_options()
      ("help,h", "print this message")
      ("random-seed,s", po::value<unsigned >(&seed)->default_value(0),  "set random seed")
      ("method,m", po::value<std::string > (&method)->default_value("bspline"), "set method, can be bspline or fbspline")
      ("order,o", po::value<unsigned > (&order)->default_value (8), "set order of interpolation")
      ("beta,b", po::value<double > (&beta)->default_value(2), "value of beta")
      ("grid-number,k", po::value<unsigned > (&ksize)->default_value(32), "number of grid on three dimesions")
      ("rcut,c", po::value<double > (&rcut)->default_value(2.5), "cut of radius")
      ("box", po::value<double > (&boxsize)->default_value(10.), "size of box on three dimesions")
      ("npart,n", po::value<unsigned > (&npart)->default_value (1000), "number of particles")
      ("config", po::value<std::string > (&config))
      ("ref-force,r", po::value<std::string > (&forceFile)->default_value ("force.ref"), "the exact force for the system")
      ("test-integer,t", po::value<int> (&testInteger)->default_value(256), "test integer")
      ("precision,p", po::value<double> (&required)->default_value(1e-3), "the target precision of optimization");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }
  if (method != std::string("bspline") && 
      method != std::string("fbspline") &&
      method != std::string("es")){
    std::cerr << "the method should be bspline or fbspline or es\ndo nothing and quit\n";
    return 0;
  }

  std::cout << "########################################\n"
	    << "# summary of parameters\n"
	    << "# method is " << method << "\n"
	    << "# order is " << order << "\n"
	    << "# beta is " << beta << "\n"
	    << "# K is " << ksize << "\n"
	    << "# rcut is " << rcut << "\n"
	    << "# box is " << boxsize << "\n"
	    << "# npart is " << npart << "\n"
	    << "# random seed is " << seed << "\n"
	    << "# configuration file is " << config << "\n"
	    << "# reference force file is " << forceFile << "\n"
	    << "# the required precision is " << required << "\n"
	    << "########################################"
	    << std::endl;
  


  // setup configuration //////////////////////////
  std::cout.precision (16);
  std::cout << std::scientific;
  double boxlx = boxsize*1;
  double boxly = boxsize*1;
  double boxlz = boxsize*1;

  std::vector<StandardParticle > particles;
  StandardParticle tmpp;
  double vaver = 10;
   
  FILE *fp = NULL;
  if (vm.count("config")){
    fp = fopen (config.c_str(), "r");
    if (fp == NULL){
      std::cerr << "cannot open file " << config << std::endl;
      return 1;
    }
    fscanf (fp, "%d", &npart);
  }
  
  ToolBox::init_genrand (seed);
  for (unsigned i = 0; i < npart; ++i){
    double theta = 2 * M_PI * ToolBox::genrand_real2();
    double phi = M_PI * ToolBox::genrand_real2();

    if (!vm.count("config")){
      tmpp.clear();
      tmpp.r()[0] = boxlx * ToolBox::genrand_real2 ();
      tmpp.r()[1] = boxly * ToolBox::genrand_real2 () ;
      tmpp.r()[2] = boxlz * ToolBox::genrand_real2 () ;
      if (ToolBox::genrand_real2() > 0.5){
	tmpp.charge() = 1;
      }
      else {
	tmpp.charge() = -1;
      }
    }
    else {
      fscanf (fp, "%lf%lf%lf%lf", &tmpp.r()[0], &tmpp.r()[1], &tmpp.r()[2], &(tmpp.charge()));
    } 

    tmpp.v()[0] = vaver*sin(phi)*cos(theta);
    tmpp.v()[1] = vaver*sin(phi)*sin(theta);
    tmpp.v()[2] = vaver*cos(phi);
    tmpp.mass() = 1.;
    tmpp.id() = i ;
    tmpp.type() = 0;

    particles.push_back (tmpp);
  }
  if (vm.count("config")){
    fscanf (fp, "%lf", &boxsize);
  }
//   boxlx = boxsize;
//   boxly = boxsize*1;
//   boxlz = boxsize*1;
  TimeConstRectangularBoxGeometry box(boxlx, boxly, boxlz);
  for (std::vector<StandardParticle >::iterator ppart = particles.begin(); 
       ppart != particles.end(); ++ppart){
    box.moveParticleToBox (*ppart);
  }
  
  if (vm.count("config")){
    std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << "\n"
	      << "$ npart in config file is " << npart << "\n"
	      << "$ boxl  in config file is " << boxsize << "\n"
	      << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << "\n";
    fclose (fp);
  }
  ////////////////////////////////////////

  std::vector<unsigned > K(3, ksize);


  double unitRequired = required / unitScale * 1.;
  EleOptimizer tc;
//   tc.init (K, rcut, box, InterpolationInfo::BSpline, order, unitRequired);
  if (method == "bspline"){
    tc.init (box, InterpolationInfo::BSpline, unitRequired); 
  }
  else if (method == "fbspline"){
    tc.init (box, InterpolationInfo::FBSpline, unitRequired); 
  }
  else {
    std::cerr << "wrong method " << std::endl;
  }  
  for (std::vector<StandardParticle >::iterator ppart = particles.begin();
       ppart != particles.end(); ++ppart){
    tc.registerParticle (*ppart);
  }
  tc.build();

  
//   double dir, loop, conv;
//   TimeCounter tc;
//   std::vector<unsigned > tmpK (3, 64);
//   tc.init(tmpK, rcut, box);
//   for (std::vector<StandardParticle >::iterator ppart = particles.begin();
//        ppart != particles.end(); ++ppart){
//     tc.registerParticle (*ppart);
//   }
//   tc.build();
//   tc.cal_time (rcut, K, order, dir, loop, conv);
//   std::cout.precision(3);
//   std::cout << "estimated dir time " << dir <<"\n"
// 	    << "estimated loop time " << loop <<"\n"
// 	    << "estimated conv time " << conv <<"\n"
//       ;
  

  
  std::vector<unsigned > newK;
  InterpolationInfo::Order neworder;
  value_type newrcut;
  value_type newbeta;
  
  tc.optiParam (newrcut, newK, neworder, newbeta, 1e-1);
  

//   std::cout << tc.findMinBeta (1e-6/unitScale) << std::endl;

//   ToolBox::IntegerDecomposer id;
//   std::vector<int > re;
//   std::cout << id.decompose (testInteger, re) << std::endl;
//   std::copy (re.begin(), re.end(), std::ostream_iterator<int >(std::cout, "\t"));
//   std::cout << std::endl;
  
  
 
  return 0;
}
