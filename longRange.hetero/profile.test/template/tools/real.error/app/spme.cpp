#include "Particle.h"
#include "ToolBox.h"
#include "ShortRangePairFinder.h"
#include "BoxGeometry.h"
#include "PairInteraction.h"
#include "VectorOperation.h"
#include "ElectrostaticInteraction.h"
#include "Integral1D.h"
#include "Polynominal.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "GroFileManager.h"

// static double unitScale = 138.935485;
static double unitScale = 1.;

int main (int argc, char * argv[])
{
  unsigned seed = 0;
  double beta;
  unsigned ksize;
  double rcut;
  double boxsize = 1;
  unsigned npart;
  std::string method;
  unsigned order;
  po::options_description desc ("Allow options");
  std::string config;
  std::string forceFile;
  std::vector<double > betaSet;
  
  desc.add_options()
      ("help,h", "print this message")
      ("method,m", po::value<std::string > (&method)->default_value("bspline"), "set method, can be bspline or fbspline")
      ("order,n", po::value<unsigned > (&order)->default_value (4), "set order of interpolation")
      ("beta,b", po::value<double > (&beta)->default_value(1.), "value of beta")
      ("grid-number,k", po::value<unsigned > (&ksize)->default_value(32), "number of grid on three dimesions")
      ("rcut,r", po::value<double > (&rcut)->default_value(2.5), "cut of radius")
      ("config,c",    po::value<std::string > (&config)->default_value ("conf.gro"), "config file")
      ("output-force,o", po::value<std::string > (&forceFile)->default_value ("force.ref"), "the exact force for the system");
  
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

  // setup configuration //////////////////////////
  std::cout.precision (16);
  std::cout << std::scientific;
  double boxlx = boxsize;
  double boxly = boxsize*1;
  double boxlz = boxsize*1;

  std::vector<StandardParticle > particles;
  StandardParticle tmpp;
  std::vector<int > resdindex, atomindex;
  std::vector<std::string > resdname, atomname;
  std::vector<std::vector<double > > posi, velo;
  std::vector<double > tmpbox;
  GroFileManager::read (config,
			resdindex, resdname, atomname, atomindex,
			posi, velo,
			tmpbox);
  npart = posi.size();

  std::vector<double > zero3 (3, 0.);
  std::vector<std::vector<double > > dir_force (npart, zero3);
  std::vector<std::vector<double > > rec_force (npart, zero3);
  
  ToolBox::init_genrand (seed);
  for (unsigned i = 0; i < npart; ++i){
    tmpp.r()[0] = posi[i][0];
    tmpp.r()[1] = posi[i][1];
    tmpp.r()[2] = posi[i][2];
    tmpp.v()[0] = 0.;
    tmpp.v()[1] = 0.;
    tmpp.v()[2] = 0.;
    tmpp.mass() = 1.;
    tmpp.id() = i ;
    tmpp.type() = 0;
    if (i < npart / 2) {
      tmpp.charge() = 1.;
    }
    else {
      tmpp.charge() = -1.;
    }
    particles.push_back (tmpp);
  }
  boxlx = tmpbox[0];
  boxly = tmpbox[1];
  boxlz = tmpbox[2];

  TimeConstRectangularBoxGeometry box(boxlx, boxly, boxlz);
  for (std::vector<StandardParticle >::iterator ppart = particles.begin(); 
       ppart != particles.end(); ++ppart){
    box.moveParticleToBox (*ppart);
  }

  std::cout.precision (6);
  std::cout << "########################################\n"
	    << "# summary of parameters\n"
	    << "# method is " << method << "\n"
	    << "# order is " << order << "\n"
	    << "# beta is " << beta << "\n"
	    << "# K is " << ksize << "\n"
	    << "# rcut is " << rcut << "\n"
	    << "# box is " << boxlx << "\n"
	    << "# npart is " << npart << "\n"
	    << "# random seed is " << seed << "\n"
	    << "# configuration file is " << config << "\n"
	    << "########################################"
	    << std::endl;
  ////////////////////////////////////////
  
  std::vector<unsigned > K(3, ksize);
  // K[0] += 2;
  // K[1] += 1;
//   double beta = 3;
  ElectrostaticInteraction ele ;
  if (method == std::string("bspline")){
    ele.init (beta, K, rcut, box, InterpolationInfo::BSpline, order);
  }
  else if (method == std::string("fbspline")){
    ele.init (beta, K, rcut, box, InterpolationInfo::FBSpline, order);
  }  
  else {
    ele.init (beta, K, rcut, box, InterpolationInfo::ES);
  }
  std::cout << "### initilized elector calculater" << std::endl;
  for (std::vector<StandardParticle >::iterator ppart = particles.begin();
       ppart != particles.end(); ppart++){
    ele.registerParticle (*ppart);
  }
  std::cout << "### registered Particle" << std::endl;

  Stopwatch tmpwatch;
  tmpwatch.start();
  ele.buildRec();
  ele.buildDir();
  tmpwatch.stop();
  std::cout.precision (3);
  std::cout << "build up time cost " << tmpwatch.user()<< std::endl;
  std::cout.precision (16);
  std::cout << "### build up calculator" << std::endl;

  double tmp0, tmp1;
  ele.errorEstimatePME (tmp0, tmp1);
  tmp0 *= unitScale;
  tmp1 *= unitScale;
  std::cout << "the solver ES dir error is " << tmp0 << std::endl;
  std::cout << "the solver ES rec error is " << tmp1 << std::endl;
  
  // clear interaction
  for (std::vector<StandardParticle >::iterator pp = particles.begin(); 
       pp != particles.end(); pp++){
    std::fill (pp->f().begin(), pp->f().end(), 0.0);
  }
  std::cout << "### force is cleared" << std::endl;
  // calculate interaction
  double potenEnergy = ele.applyInteractionCalPotential (dir_force, rec_force);
  ele.reset_watch();
  Stopwatch st;
  st.start();
  // for (int i = 0; i < 1; ++i){
  //   ele.applyInteraction();
  // }
  st.stop();
  double dirtime, rectime, looptime, convtime;
  ele.get_time (dirtime, rectime, looptime, convtime);
  std::cout.precision (3);
  double scale = 1./1.;
  std::cout << "time rec " << rectime*scale << std::endl;
  std::cout << "time rec loop " << looptime*scale << std::endl;
  std::cout << "time rec conv " << convtime*scale << std::endl;
  std::cout << "time dir " << dirtime*scale << std::endl;
  std::cout << "time total " << (st.user())*scale << std::endl;
  std::cout.precision (16);
  std::cout << "potential energy is " << potenEnergy << std::endl;
  std::cout << "### interaction is applied" << std::endl;
  
//   ele.calPotential (
//       allpair.begin(), allpair.end());
//   refEle.calPotential (
//       refAllpair.begin(), refAllpair.end());

//   std::cout << erfc(beta*1.41) / rcut << std::endl;
  
//   std::cerr.precision(10);
//   std::cerr <<potenEnergy << std::endl;
//   std::cerr << particles[0].f()[0] << '\t'
// 	    << particles[0].f()[1] << '\t'
// 	    << particles[0].f()[2] << '\n';

  double errorRec;
  ele.errorEstimateInterpolRec (errorRec);
  errorRec *= unitScale;
  std::cout << "the estimated error of rec part SPME " << errorRec << std::endl;

  FILE * fpmy = fopen (forceFile.c_str(), "w");
  if (fpmy == NULL){
    std::cerr << "cannot open file " << forceFile << std::endl;
    return 1;
  }
  for (unsigned i = 0; i < particles.size(); ++ i){
    fprintf (fpmy, "%.16e %.16e %.16e \t  %.16e %.16e %.16e \t  %.16e %.16e %.16e\n", 
	     dir_force[i][0],
	     dir_force[i][1],
	     dir_force[i][2],
	     unitScale * particles[i].f()[0],
	     unitScale * particles[i].f()[1],
	     unitScale * particles[i].f()[2],
	     dir_force[i][0] + rec_force[i][0],
	     dir_force[i][1] + rec_force[i][1],
	     dir_force[i][2] + rec_force[i][2]
	);
    // fprintf (fpmy, "%.16e\t%.16e\t%.16e\n", 
    // 	     unitScale * particles[i].f()[0],
    // 	     unitScale * particles[i].f()[1],
    // 	     unitScale * particles[i].f()[2]);
  }
  fclose (fpmy);

  return 0;

}




//   std::vector<PositionParticle *> tmpparts;
//   tmpparts.push_back(&particles[0]);
//   tmpparts.push_back(&particles[1]);
//   tmpparts.push_back(&particles[2]);
  
//   allpair.unregisterParticle (tmpparts);
//   allpair.rebuild();
//   allpair.registerParticle (tmpparts);
//   allpair.rebuild();
//   for (std::vector<PositionParticle * >::iterator itList = allpair.begin();
//        itList != allpair.end(); ){
//     int a = (*(itList++))->id();
//     int b = (*(itList++))->id();
//     if ( a == 0 || b == 0)
//       std::cout << a << '\t'
// 		<< b << '\n';
//   }




//   for (int nt = 0; nt < nstep; nt ++){
//     for (std::vector<StandardParticle >::iterator pp = particles.begin(); pp != particles.end(); pp++){
//       VectorOperation::add (pp->v(), 0.5/pp->mass()*dt, pp->f());
//       VectorOperation::add (pp->r(), dt, pp->v());
//     }
//     // move particles back to box
//     box.moveParticleToBox (time, particles);
      
//     // clear interaction
//     for (std::vector<StandardParticle >::iterator pp = particles.begin(); pp != particles.end(); pp++){
//       std::fill (pp->f().begin(), pp->f().end(), 0.0);
//     }
//     // calculate interaction
//     potenEnergy = 0;
//     int count = 0;
//     allpair.reset();
//     while (allpair.getPair (p0, p1)){
//       count ++;
//       potenEnergy +=
// 	  lj.applyInteractionCalPotential (time, 
// 					   * static_cast<StandardParticle * > (p0),
// 					   * static_cast<StandardParticle * > (p1), 
// 					   static_cast<BoxGeometry &> (box));
//     }
//     // step forward
//     for (std::vector<StandardParticle >::iterator pp = particles.begin(); pp != particles.end(); pp++){
//       VectorOperation::add (pp->v(), 0.5/pp->mass()*dt, pp->f());
//     }
//     // static kinetic energy
//     kineticEnergy = 0;
//     for (std::vector<StandardParticle >::iterator pp = particles.begin(); pp != particles.end(); pp++){
//       kineticEnergy += 0.5 * VectorOperation::dot (pp->v(), pp->v());
//     }

//     time += dt;
//     // print
//     if (nt % 100 == 0)
//       std::cout << time << '\t'
// 		<< kineticEnergy << '\t'
// 		<< potenEnergy << '\t'
// 		<< kineticEnergy + potenEnergy <<'\n';
  
//   }
  
