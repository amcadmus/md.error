#include "Particle.h"
#include "ToolBox.h"
#include "ShortRangePairFinder.h"
#include "BoxGeometry.h"
#include "PairInteraction.h"
#include "VectorOperation.h"
#include "ElectrostaticInteraction.h"
#include "Integral1D.h"
#include "Polynominal.h"

// static double unitScale = 138.935485;
static double unitScale = 1.;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main (int argc, char * argv[])
{
  unsigned seed;
  double refRcut;
  double boxsize;
  unsigned npart;
  po::options_description desc ("Allow options");
  std::string config;
  std::string forceFile;
  double precision;
  double definedbeta;
  unsigned definedk;
  
  desc.add_options()
      ("help,h", "print this message")
      ("random-seed,s", po::value<unsigned >(&seed)->default_value(0),  "set random seed")
      ("box", po::value<double > (&boxsize)->default_value(10.), "size of box on three dimesions")
      ("npart,n", po::value<unsigned > (&npart)->default_value (1000), "number of particles")
      ("config,c", po::value<std::string > (&config))
      ("beta,b", po::value<double > (&definedbeta)->default_value(2), "value of beta")
      ("grid-number,k", po::value<unsigned > (&definedk)->default_value(32), "number of grid on three dimesions")
      ("precision,p", po::value<double > (&precision), "required precision of the reference force. if this option is used, then beta and grid-number are ignored")
      ("ref-force,r", po::value<std::string > (&forceFile)->default_value ("force.ref"), "the output exact force for the system");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  if (vm.count("precision")){
    std::cout << "##############################\n"
	      << "# summary of parameters\n"
	      << "# box is " << boxsize << "\n"
	      << "# npart is " << npart << "\n"
	      << "# random seed is " << seed << "\n"
	      << "# configuration file is " << config << "\n"
	      << "# required precision is " << precision << "\n"
	      << "# output reference force file is " << forceFile << "\n"
	      << "##############################\n"
	      << std::endl;
  }
  else{
    std::cout << "##############################\n"
	      << "# summary of parameters\n"
	      << "# box is " << boxsize << "\n"
	      << "# npart is " << npart << "\n"
	      << "# random seed is " << seed << "\n"
	      << "# configuration file is " << config << "\n"
	      << "# beta is " << definedbeta << "\n"
	      << "# K is " << definedk << "\n"
	      << "# output reference force file is " << forceFile << "\n"
	      << "##############################\n"
	      << std::endl;
  }
  

  // setup configuration //////////////////////////
  std::cout.precision (16);
  std::cout << std::scientific;
  double boxlx = boxsize;
  double boxly = boxsize*1;
  double boxlz = boxsize*1;

  std::vector<StandardParticle > refParticles;
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

    refParticles.push_back (tmpp);
  }
  if (vm.count("config")){
    fscanf (fp, "%lf", &boxsize);
  }
  boxlx = boxsize;
  boxly = boxsize*1;
  boxlz = boxsize*1;
  std::cout << boxlx << std::endl;

  refRcut = 0.9 * (boxsize / 2.);
  if (refRcut > 5){
    refRcut = 5;
  }  


  TimeConstRectangularBoxGeometry box(boxlx, boxly, boxlz);
  for (std::vector<StandardParticle >::iterator ppart = refParticles.begin(); 
       ppart != refParticles.end(); ++ppart){
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

  
  double c2 = 0;
  int Npart = refParticles.size();
  for (std::vector<StandardParticle >::iterator ppart = refParticles.begin(); ppart != refParticles.end(); ppart++){
    c2 += (ppart)->charge() * (ppart)->charge();
  }
  std::vector<unsigned > refK (3, 10);
//   double beta = 3;
  ElectrostaticInteraction refEle;
  double refbeta ;
  double tmp0;
  double tmp1;

//   std::vector<double > tmp;
//   box.getBoxSize(tmp);
//   std::cout << tmp[0] << std::endl;
//   std::cout << tmp[1] << std::endl;
//   std::cout << tmp[2] << std::endl;
  
  if (vm.count("precision")){
    refbeta = 0.5;
//     std::cout << refbeta << '\t'
// 	      << rdfK << '\t'
// 	      << refRcut << '\t'
// 	      << std::endl;
    refEle.init (refbeta, refK, refRcut, box, InterpolationInfo::ES);
    std::cout << "### initilized elector calculater" << std::endl;
    refEle.errorEstimatePME (c2, Npart, tmp0, tmp1);
    tmp0 *= unitScale;
    tmp1 *= unitScale;
    while (tmp0 > 0.5 * precision) {
      if (refbeta > 2.5) break;
      refbeta += 0.2;
      refEle.reinit (refbeta, refK, refRcut, box, InterpolationInfo::ES);
      refEle.errorEstimatePME (c2, Npart, tmp0, tmp1);
      tmp0 *= unitScale;
      tmp1 *= unitScale;
//     std::cout << "**** Warning *** the rcut of real space is too small \n" ;
    }
    while (tmp1 > 0.5 * precision) {
      if (refK[0] >= 200) break;
      refK[0] += 10;
      refK[1] += 10;
      refK[2] += 10;
      refEle.reinit (refbeta, refK, refRcut, box, InterpolationInfo::ES);
      refEle.errorEstimatePME (c2, Npart, tmp0, tmp1);
      tmp0 *= unitScale;
      tmp1 *= unitScale;
    }
  }
  else {
    refbeta = definedbeta;
    refK[0] = definedk;
    refK[1] = definedk;
    refK[2] = definedk;
    refEle.init (refbeta, refK, refRcut, box, InterpolationInfo::ES);
    std::cout << "### initilized elector calculater" << std::endl;
    refEle.errorEstimatePME (c2, Npart, tmp0, tmp1);
    tmp0 *= unitScale;
    tmp1 *= unitScale;
  }
  std::cout << "ref beta is " << refbeta << std::endl;
  std::cout << "ref K0 is " << refK[0] << std::endl;
  std::cout << "size of K is " << refK[0] << std::endl;
  std::cout << "the ref ES dir error is " << tmp0 << std::endl;
  std::cout << "the ref ES rec error is " << tmp1 << std::endl;

  for (std::vector<StandardParticle >::iterator ppart = refParticles.begin(); ppart != refParticles.end(); ppart++){
    refEle.registerParticle(*ppart);
  }
  std::cout << "### registed particles\n" ;
  refEle.build();
  std::cout << "### build all\n";


  // clear interaction
  for (std::vector<StandardParticle >::iterator pp = refParticles.begin(); 
       pp != refParticles.end(); pp++){
    std::fill (pp->f().begin(), pp->f().end(), 0.0);
  }
  std::cout << "### force is cleared" << std::endl;
  // calculate interaction
  double refPotenEnergy = 0;
  std::cout << "### list checked" << std::endl;
  refPotenEnergy += refEle.applyInteractionCalPotential ();
  std::cout << "### interaction is applied" << std::endl;
  
  std::cout << "reference energy is " << refPotenEnergy << std::endl;

  // print reference force
  fp = fopen (forceFile.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << forceFile << std::endl;
    return 1;
  }
  for (unsigned i = 0; i < refParticles.size(); ++ i){
    fprintf (fp, "%.16e\t%.16e\t%.16e\n", 
	     unitScale * refParticles[i].f()[0],
	     unitScale * refParticles[i].f()[1],
	     unitScale * refParticles[i].f()[2]);
  }
  

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
  
