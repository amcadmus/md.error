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

// static double unitScale = 138.935485;
static double unitScale = 1.;

int main (int argc, char * argv[])
{
  unsigned seed;
  double beta;
  unsigned ksize;
  double rcut;
  double refRcut;
  double boxsize;
  unsigned npart;
  std::string method;
  unsigned order;
  po::options_description desc ("Allow options");
  std::string config;
  std::string forceFile;
  std::vector<double > betaSet;
  
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
      ("beta-set", po::value<std::vector<double > > (), "a set of values of beta, --beta will be ignored");
  
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
	    << "########################################"
	    << std::endl;
  refRcut = boxsize / 2.;
  
  // load reference force /////////////////////////
  FILE * fp = fopen(forceFile.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << forceFile << std::endl;
    return 1;
  }
  std::vector<std::vector<double > > refForce;
  std::vector<double > tmpF (3);
  while (fscanf(fp, "%lf%lf%lf", &tmpF[0], &tmpF[1], &tmpF[2]) != EOF){
    refForce.push_back(tmpF);
  }
  fclose (fp);
  //////////////////////////////////////////////////



  // setup configuration //////////////////////////
  std::cout.precision (16);
  std::cout << std::scientific;
  double boxlx = boxsize;
  double boxly = boxsize*1;
  double boxlz = boxsize*1;

  std::vector<StandardParticle > particles;
  StandardParticle tmpp;
  double vaver = 10;
  double potenEnergy = 0;
   
  fp = NULL;
  if (vm.count("config")){
    fp = fopen (config.c_str(), "r");
    if (fp == NULL){
      std::cerr << "cannot open file " << config << std::endl;
      return 1;
    }
    fscanf (fp, "%d", &npart);
  }
  if (npart != refForce.size()){
    std::cerr <<"the size of the system should be the same as the number of reference force." << std::endl;
    return 1;
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
  boxlx = boxsize;
  boxly = boxsize*1;
  boxlz = boxsize*1;
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
  
  if (!vm.count("beta-set")){
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
    for (std::vector<StandardParticle >::iterator ppart = particles.begin(); ppart != particles.end(); ppart++){
      ele.registerParticle (*ppart);
    }
    std::cout << "### registered Particle" << std::endl;
//     double a, b, c;
//     ele.errorEstimatePME (a, b);
//     ele.errorEstimateInterpolRec (c);
//     std::cout << "pre build error estimates\n"
// 	      << a*unitScale << "\n"
// 	      << b*unitScale << "\n"
// 	      << c*unitScale << "\n";
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
    potenEnergy = ele.applyInteractionCalPotential ();
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

    std::vector<StandardParticle >::iterator ppart = particles.begin();
    std::vector<std::vector<double > >::iterator prefForce = refForce.begin();
    double sum = 0;
    for (; ppart != particles.end(); ++ppart, ++prefForce){
      std::vector<double > diff (ppart->f());
      VectorOperation::scale (diff, unitScale);
      VectorOperation::add (diff, -1, *prefForce);
      sum += VectorOperation::dot (diff, diff);
    }
    sum /= particles.size();
    std::cout << "average diff sqruare is " << (sum) 
	      << " and its square root is " << sqrt (sum)
	      << std::endl;
    double errorRec;
    ele.errorEstimateInterpolRec (errorRec);
    errorRec *= unitScale;
    std::cout << "the estimated error of rec part SPME " << errorRec << std::endl;
  }
  else{
    std::vector<unsigned > K(3, ksize);
    ElectrostaticInteraction ele ;
    bool firstTime = true;
    betaSet = vm["beta-set"].as< std::vector<double> >();
    
    std::cout << betaSet.size() << std::endl;
    for (std::vector<double >::iterator pbeta = betaSet.begin();
	 pbeta != betaSet.end(); ++pbeta){
      std::cout << "### using beta value " << *pbeta << std::endl;
      if (firstTime){
	if (method == std::string("bspline")){
	  ele.init (*pbeta, K, rcut, box, InterpolationInfo::BSpline, order);
	}
	else if (method == std::string("fbspline")){
	  ele.init (*pbeta, K, rcut, box, InterpolationInfo::FBSpline, order);
	}  
	else {
	  ele.init (*pbeta, K, rcut, box, InterpolationInfo::ES);
	}
      }
      else {
	if (method == std::string("bspline")){
	  ele.reinit (*pbeta, K, rcut, box, InterpolationInfo::BSpline, order);
	}
	else if (method == std::string("fbspline")){
	  ele.reinit (*pbeta, K, rcut, box, InterpolationInfo::FBSpline, order);
	}  
	else {
	  ele.reinit (*pbeta, K, rcut, box, InterpolationInfo::ES);
	}
      }
      firstTime = false;
      
      for (std::vector<StandardParticle >::iterator ppart = particles.begin(); ppart != particles.end(); ppart++){
	ele.registerParticle (*ppart);
      }
      
      ele.buildRec();
      ele.buildDir();
    
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
      Stopwatch st;
      ele.reset_watch();
      st.start();
      for (int i = 0; i < 1; ++i){
	ele.applyInteraction();
      }
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
//       std::cout << "potential energy is " << potenEnergy << std::endl;
//       std::cout << "### interaction is applied" << std::endl;

      std::vector<StandardParticle >::iterator ppart = particles.begin();
      std::vector<std::vector<double > >::iterator prefForce = refForce.begin();
      double sum = 0;
      for (; ppart != particles.end(); ++ppart, ++prefForce){
	std::vector<double > diff (ppart->f());
	VectorOperation::scale (diff, unitScale);
	VectorOperation::add (diff, -1, *prefForce);
	sum += VectorOperation::dot (diff, diff);
      }
      sum /= particles.size();
      std::cout << "average diff sqruare is " << (sum) 
		<< " and its square root is " << sqrt (sum)
		<< std::endl;
      double errorRec;
      ele.errorEstimateInterpolRec (errorRec);
      errorRec *= unitScale;
      std::cout << "the estimated error of rec part SPME " << errorRec << std::endl;
    }
  }  

  FILE * fpmy = fopen ("myForce", "w");
  for (unsigned i = 0; i < particles.size(); ++ i){
    fprintf (fpmy, "%.16e\t%.16e\t%.16e\n", 
	     unitScale * particles[i].f()[0],
	     unitScale * particles[i].f()[1],
	     unitScale * particles[i].f()[2]);
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
  
