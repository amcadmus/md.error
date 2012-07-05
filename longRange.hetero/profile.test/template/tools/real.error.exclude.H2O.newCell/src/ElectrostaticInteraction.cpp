#include "ElectrostaticInteraction.h"
#include "Polynominal.h"
#include "VectorOperation.h"
#include "Integral3D.h"
#include <numeric>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include "CellList.h"

// void print3D (const std::string & filename, 
// 	      const int & K0, const int & K1, const int & K2,
// 	      fftw_complex * v)
// {
//   FILE * fp = fopen (filename.c_str(), "w");
//   for (int i = 0; i < K0; i ++){
//     for (int j = 0; j < K1; j ++){
//       for (int k = 0; k < K2; k ++){
// 	fprintf (fp, "%.6e ", v[k+K2*(j+K1*i)][0]);
//       }
//       fprintf (fp, "\n");
//     }
//     fprintf (fp, "\n");
//   }
//   fclose (fp);
// }
struct Compare : public std::binary_function<StandardParticle *, StandardParticle *, bool> {
  bool operator()(StandardParticle * x, StandardParticle * y){ 
    return x->id() < y->id(); 
  }
};

ElectrostaticInteraction::ElectrostaticInteraction()
{  
  dir_time = 0;
  rec_time = 0;
}

void ElectrostaticInteraction::init (const value_type & beta_,
				     const std::vector<unsigned > & K_,
				     const value_type & rcut_, 
				     const TimeConstRectangularBoxGeometry & box_,
				     const InterpolationInfo::Method & method_,
				     const InterpolationInfo::Order  & order_)
{
  beta = beta_;
  K = K_;
  rcut = rcut_;
  box = &box_;
  method = method_;
  order = order_;
  box->getBaseVector (baseVector);
  calV();
  calAStar();
  
  dir = new ElectrostaticInteraction_dir (rcut, beta);
  pairFinder = new NeighborListPairFinder (rcut, 0.1 * rcut, *box);
  
  if (method == InterpolationInfo::ES){
    rec = new ElectrostaticInteraction_rec_ES ();
    static_cast<ElectrostaticInteraction_rec_ES * >(rec) -> init (baseVector, K, beta);
  }
  else if (method == InterpolationInfo::BSpline){
    rec = new ElectrostaticInteraction_rec_BSpline ();
    static_cast<ElectrostaticInteraction_rec_BSpline * >(rec) -> init (baseVector, K, beta, order);
  }
  else if (method == InterpolationInfo::FBSpline){
    rec = new ElectrostaticInteraction_rec_FBSpline ();
    static_cast<ElectrostaticInteraction_rec_FBSpline * >(rec) -> init (baseVector, K, beta, order);
  }
  else {
    std::cerr << "no such method" << method 
	      << "\nuse Bspline" << std::endl;
    rec = new ElectrostaticInteraction_rec_BSpline ();
  }
}

void ElectrostaticInteraction::buildRec ()
{  
  if (method == InterpolationInfo::ES){
  }
  else if (method == InterpolationInfo::BSpline){
    static_cast<ElectrostaticInteraction_rec_BSpline * >(rec) -> build();
  }
  else if (method == InterpolationInfo::FBSpline){
    static_cast<ElectrostaticInteraction_rec_FBSpline * >(rec) -> build();
  }
  else {
    static_cast<ElectrostaticInteraction_rec_BSpline * >(rec) -> build();
  }
}

void ElectrostaticInteraction::registerParticle (StandardParticle & p)
{
  partPool.push_back (&p);
  // pairFinder->registerParticle (p);
}

void ElectrostaticInteraction::buildDir ()
{
  // pairFinder->build();
}

void ElectrostaticInteraction::build ()
{
  buildRec();
  buildDir();
}

void ElectrostaticInteraction::unregisterParticle (StandardParticle & p)
{
  for (std::vector<StandardParticle *>::iterator i = partPool.begin();
       i != partPool.end(); ++ i){
    if ((*i)->id() == p.id()){
      partPool.erase(i);
    }
  }
  // pairFinder->unregisterParticle (p);
}

void ElectrostaticInteraction::rebuildDir ()
{
  // pairFinder->rebuild();
}

void ElectrostaticInteraction::clear () 
{
  partPool.clear();
  delete dir;
  delete rec;
  delete pairFinder;
}


void ElectrostaticInteraction::applyInteraction ()
{
  watch.start();
  // pairFinder->checkList();
  // for (std::vector<PositionParticle * >::const_iterator ppart = pairFinder->begin();
  //      ppart != pairFinder->end(); ){
  //   dir->applyInteraction (
  // 	* static_cast<StandardParticle * >(*(ppart++)), 
  // 	* static_cast<StandardParticle * >(*(ppart++)),
  // 	* box);
  // }
  watch.stop();
  dir_time += watch.user();
  
  watch.start();
  rec->applyInteraction (partPool.begin(), partPool.end());
  watch.stop();
  rec_time += watch.user();
}

value_type ElectrostaticInteraction::applyInteractionCalPotential ()
{
  double valueRec = 0;
  double valueDir = 0;
  double valueCor = 0;

  watch.start();
  // pairFinder->checkList();
  // for (std::vector<PositionParticle * >::const_iterator ppart = pairFinder->begin();
  //      ppart != pairFinder->end(); ){
  //   valueDir += dir->applyInteractionCalPotential (
  // 	* static_cast<StandardParticle * >(*(ppart++)), 
  // 	* static_cast<StandardParticle * >(*(ppart++)),
  // 	* box);
  // }
  watch.stop();
  dir_time += watch.user();  
  
  watch.start();
  // for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
  //      ppart != partPool.end(); ppart ++){
  //   valueCor -= beta / sqrt(M_PI) * (*ppart)->charge() * (*ppart)->charge();
  // }
  watch.stop();
  corr_time += watch.user();
  
  watch.start();
  valueRec = rec->applyInteractionCalPotential (partPool.begin(), partPool.end());
  watch.stop();
  rec_time += watch.user();

  return valueCor + valueDir + valueRec;
}

value_type ElectrostaticInteraction::
applyInteractionCalPotential (std::vector<std::vector<double > > & dir_force,
			      std::vector<std::vector<double > > & rec_force)
{
  double valueRec = 0;
  double valueDir = 0;
  double valueCor = 0;

  // watch.start();
  // pairFinder->checkList();
  // for (std::vector<PositionParticle * >::const_iterator ppart = pairFinder->begin();
  //      ppart != pairFinder->end(); ){
  //   valueDir += dir->applyInteractionCalPotential (
  // 	* static_cast<StandardParticle * >(*(ppart++)), 
  // 	* static_cast<StandardParticle * >(*(ppart++)),
  // 	* box);
  // }
  // watch.stop();
  // dir_time += watch.user();
  
  std::vector<double > tmpmybox;
  box->getBoxSize (tmpmybox);
  VectorType mybox;
  mybox.x = tmpmybox[0];
  mybox.y = tmpmybox[1];
  mybox.z = tmpmybox[2];
  std::vector<std::vector<double > > coords;
  std::vector<double > charges;
  for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
       ppart != partPool.end(); ppart ++){
    coords.push_back ((*ppart)->r());
    for (unsigned dd = 0; dd < 3; ++dd){
      if (coords.back()[dd] < 0.) {
	coords.back()[dd] += tmpmybox[dd];
      }
      else if (coords.back()[dd] >= tmpmybox[dd]){
	coords.back()[dd] -= tmpmybox[dd];
      }
    }
    charges.push_back ((*ppart)->charge());
  }
  
  CellList clist (coords.size(), mybox, rcut);
  printf ("### build cell list\n");
  clist.rebuild (coords);
  printf ("### cal dir interaction\n");
  for (unsigned ii = 0; ii < clist.getTotalNumCell(); ++ii){
    printf ("### %d-th cell, %d in total   \r", ii+1, clist.getTotalNumCell());
    fflush (stdout);      
    std::vector<unsigned > ilist = clist.getList()[ii];
    std::vector<unsigned > neighborCellIndex = clist.neighboringCellIndex (ii);
    // if (ii == 0){
    //   for (unsigned jj = 0; jj < neighborCellIndex.size(); ++jj){
    // 	printf ("neighbor %d cell: %d\n", jj, neighborCellIndex[jj]);
    //   }
    // }
    for (unsigned jj = 0; jj < neighborCellIndex.size(); ++jj){
      if (ii > neighborCellIndex[jj]) continue;
      std::vector<unsigned > jlist = clist.getList()[neighborCellIndex[jj]];
      for (unsigned iIdx = 0; iIdx < ilist.size(); ++iIdx){
	unsigned iatom = ilist[iIdx];
	for (unsigned jIdx = 0; jIdx < jlist.size(); ++jIdx){
	  unsigned jatom = jlist[jIdx];
	  if (ii == neighborCellIndex[jj] && iatom >= jatom) continue;
	  // if (iatom == jatom) continue;
	  std::vector<double > f0;
	  dir->applyInteractionCalPotential (coords[iatom], coords[jatom], charges[iatom], charges[jatom], f0, *box);
	  for (unsigned dd = 0; dd < 3; ++dd){
	    dir_force[iatom][dd] += f0[dd];
	    dir_force[jatom][dd] -= f0[dd];
	  }
	}
      }
    }
  }
  printf ("\n");

  // unsigned count = 0;
  // for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
  //      ppart != partPool.end(); ppart ++){
  //   dir_force[count ++] = (*ppart)->f();
  //   (*ppart)->f()[0] = 0.;
  //   (*ppart)->f()[1] = 0.;
  //   (*ppart)->f()[2] = 0.;
  // }

  for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
       ppart != partPool.end(); ppart ++){
    (*ppart)->f()[0] = 0.;
    (*ppart)->f()[1] = 0.;
    (*ppart)->f()[2] = 0.;
  }
  
  watch.start();
  for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
       ppart != partPool.end(); ppart ++){
    valueCor -= beta / sqrt(M_PI) * (*ppart)->charge() * (*ppart)->charge();
  }
  watch.stop();
  corr_time += watch.user();
  
  watch.start();
  valueRec = rec->applyInteractionCalPotential (partPool.begin(), partPool.end());
  watch.stop();
  rec_time += watch.user();

  unsigned count = 0;
  for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
       ppart != partPool.end(); ppart ++){
    rec_force[count ++] = (*ppart)->f();
  }
  
  return valueCor + valueDir + valueRec;
}


value_type ElectrostaticInteraction::calPotential ()
{
  double valueRec = 0;
  double valueDir = 0;
  double valueCor = 0;

  watch.start();
  valueRec = rec->calPotential (partPool.begin(), partPool.end());
  watch.stop();
  rec_time += watch.user();
  
  watch.start();
  // pairFinder->checkList();
  // for (std::vector<PositionParticle * >::const_iterator ppart = pairFinder->begin();
  //      ppart != pairFinder->end(); ){
  //   double tmp = dir->calPotential (
  // 	* static_cast<StandardParticle * >(*(ppart++)), 
  // 	* static_cast<StandardParticle * >(*(ppart++)),
  // 	* box);
  //   valueDir += tmp;
  // }
  watch.stop();
  dir_time += watch.user();
  
  watch.start();
  for (std::vector<StandardParticle * >::const_iterator ppart = partPool.begin();
       ppart != partPool.end(); ppart ++){
    valueCor -= beta / sqrt(M_PI) * (*ppart)->charge() * (*ppart)->charge();
  }
  watch.stop();
  corr_time += watch.user();

  return valueCor + valueDir + valueRec;
}

void ElectrostaticInteraction::errorEstimateInterpolRec (value_type & error)
{
  error = -1;
  if (method == InterpolationInfo::BSpline){
    error = dynamic_cast <ElectrostaticInteraction_rec_BSpline * > 
	(rec)->errorEstimate (partPool.begin(), partPool.end());
  }
  else if (method == InterpolationInfo::FBSpline){
    error = dynamic_cast <ElectrostaticInteraction_rec_FBSpline * > 
	(rec)->errorEstimate (partPool.begin(), partPool.end());
  } 
//   std::cout << "the estimate of the dir error is " 
// 	    << dynamic_cast <ElectrostaticInteraction_rec_BSpline * > (rec)->errorEstimateDir (rcut) << std::endl;
}


void ElectrostaticInteraction::errorEstimateInterpolRec (
    const value_type & c2, const value_type & N,
    value_type & error)
{
  error = -1;
  if (method == InterpolationInfo::BSpline){
    error = dynamic_cast <ElectrostaticInteraction_rec_BSpline * > 
	(rec)->errorEstimate (c2, N);
  }
  else if (method == InterpolationInfo::FBSpline){
    error = dynamic_cast <ElectrostaticInteraction_rec_FBSpline * >
	(rec)->errorEstimate (c2, N);
  } 
//   std::cout << "the estimate of the dir error is " 
// 	    << dynamic_cast <ElectrostaticInteraction_rec_BSpline * > (rec)->errorEstimateDir (rcut) << std::endl;
}


double global_beta;
std::vector<std::vector<double > > global_vecAStar;

typedef double (*F3) (double, double, double );

double recErrorFunc (double m0, double m1, double m2)
{
  std::vector<double > m(3, 0);
  m[0] = m0 * global_vecAStar[0][0] + m1 * global_vecAStar[1][0] + m2 * global_vecAStar[2][0];
  m[1] = m0 * global_vecAStar[0][1] + m1 * global_vecAStar[1][1] + m2 * global_vecAStar[2][1];
  m[2] = m0 * global_vecAStar[0][2] + m1 * global_vecAStar[1][2] + m2 * global_vecAStar[2][2];
  double mm = VectorOperation::dot (m, m);
  return exp (- 2*M_PI*M_PI/global_beta/global_beta*mm) / mm;
}  

ElectrostaticInteraction::integralTarget::integralTarget ( ElectrostaticInteraction & ele)
    : beta (ele.beta), vecAStar (ele.vecAStar)
{
//   std::cout << "integral target " << beta << '\n';
//   std::copy (vecAStar[2].begin(), vecAStar[2].end(), std::ostream_iterator<double >(std::cout, "\n"));
}

double ElectrostaticInteraction::integralTarget::operator ()
    (double x0, double x1, double x2) const
{
  double m0 = 1./x0;
  double m1 = 1./x1;
  double m2 = 1./x2;
  std::vector<double > m(3);
  
  m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
  m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
  m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
  double mm = VectorOperation::dot (m, m);
  return m0*m0 * m1*m1 * m2*m2 * exp (-2 * M_PI * M_PI / beta / beta * mm) / mm;
  
}


void ElectrostaticInteraction::errorEstimatePME (
    const value_type & c2, const value_type & N, 
    value_type & errorDir, value_type & errorRec)
{
  errorDir = 2 * sqrt(c2/double(N))
      * sqrt (c2 / rcut / V) 
      * exp (- beta*beta*rcut*rcut);

// version 3
  double sum = 0;
  std::vector<double > m(3);
  int count (0);
  for (int m0 = -int(K[0]); m0 <= int(K[0]); ++m0){
    for (int m1 = -int(K[1]); m1 <= int(K[1]); ++m1){
      for (int m2 = -int(K[2]); m2 <= int(K[2]); ++m2){
	if (m0 >= - int(K[0])/2 && m0 <= int(K[0])/2 &&
	    m1 >= - int(K[1])/2 && m1 <= int(K[1])/2 &&
	    m2 >= - int(K[2])/2 && m2 <= int(K[2])/2) continue;
	m[0] = m0 * vecAStar[0][0] + m1 * vecAStar[1][0] + m2 * vecAStar[2][0];
	m[1] = m0 * vecAStar[0][1] + m1 * vecAStar[1][1] + m2 * vecAStar[2][1];
	m[2] = m0 * vecAStar[0][2] + m1 * vecAStar[1][2] + m2 * vecAStar[2][2];
	double mm = VectorOperation::dot (m, m);
	sum += exp (-2 * M_PI * M_PI / beta / beta * mm) / mm;
	count ++;
      }
    }
  }
//   std::cout << "k0 is " << K[0] << std::endl;
//   std::cout << "coutn is " << count << std::endl;
  errorRec =  sqrt(sum / double(N)) * 2 * c2 / V;
  

// // version 2
//   ToolBox::Integral3D<integralTarget, double > inte;
//   integralTarget func (*this);
//   double tmp= inte.cal_int (ToolBox::Integral3DInfo::Gauss27, func,
// 			    -2. / double(K[0]-4),  2. / double(K[0]-4),
// 			    -2. / double(K[1]-4),  2. / double(K[1]-4),
// 			    -2. / double(K[2]-4),  2. / double(K[2]-4),
// 			    1e-6, 20, 20, 20);
//   errorRec = sqrt(tmp / double(N)) * 2 * c2 / V;
//   std::cout << 2. / double(K[0]-4) << std::endl;
//   std::cout << tmp << std::endl;


// version 1
//   global_beta = beta;
//   global_vecAStar = vecAStar;
  
//   errorRec = sqrt(c2/double(N)) * beta / pow(V, 1./3.) / M_PI
//       * sqrt (8 * c2 / double(K[0]/2.)) 
//       * exp (-(M_PI*K[0]/2./beta/pow(V, 1./3.))*(M_PI*K[0]/2./beta/pow(V, 1./3.)));
}


void ElectrostaticInteraction::errorEstimatePME (value_type & errorDir,
						 value_type & errorRec)
{  
  std::vector<double > qs;
  for (std::vector<StandardParticle * >::iterator pp = partPool.begin();
       pp != partPool.end(); ++pp){
    qs.push_back ((*pp)->charge());
  }
  std::sort (qs.begin(), qs.end());
  // double qmax = qs[qs.size()-1];
  double Q = 0;
  for (std::vector<double >::iterator i = qs.begin();
       i != qs.end(); ++ i){
    Q += *i * *i;
  }
  errorEstimatePME (Q, partPool.size(), errorDir, errorRec);
}

void ElectrostaticInteraction::calV()
{
  V = baseVector[0][0] * (baseVector[1][1]*baseVector[2][2] - baseVector[2][1]*baseVector[1][2]) - 
      baseVector[0][1] * (baseVector[1][0]*baseVector[2][2] - baseVector[2][0]*baseVector[1][2]) +
      baseVector[0][2] * (baseVector[1][0]*baseVector[2][1] - baseVector[2][0]*baseVector[1][1]);
}
  
void ElectrostaticInteraction::calAStar ()
{
  vecAStar.resize (3);
  vecAStar[0].resize (3);
  vecAStar[1].resize (3);
  vecAStar[2].resize (3);
  vecAStar[0][0] =( baseVector[1][1]*baseVector[2][2] - baseVector[2][1]*baseVector[1][2]) / V;
  vecAStar[1][1] =( baseVector[0][0]*baseVector[2][2] - baseVector[2][0]*baseVector[0][2]) / V;
  vecAStar[2][2] =( baseVector[0][0]*baseVector[1][1] - baseVector[1][0]*baseVector[0][1]) / V;
  vecAStar[1][0] =(-baseVector[1][0]*baseVector[2][2] + baseVector[2][0]*baseVector[1][2]) / V;
  vecAStar[2][0] =( baseVector[1][0]*baseVector[2][1] - baseVector[2][0]*baseVector[1][1]) / V;
  vecAStar[0][1] =(-baseVector[0][1]*baseVector[2][2] + baseVector[2][1]*baseVector[0][2]) / V;
  vecAStar[2][1] =(-baseVector[0][0]*baseVector[2][1] + baseVector[2][0]*baseVector[0][1]) / V;
  vecAStar[0][2] =( baseVector[0][1]*baseVector[1][2] - baseVector[1][1]*baseVector[0][2]) / V;
  vecAStar[1][2] =(-baseVector[0][0]*baseVector[1][2] + baseVector[1][0]*baseVector[0][2]) / V;
}


// void BestBetaFinder::init (const std::vector<unsigned > & K_,
// 			   const value_type & rcut_, 
// 			   const TimeConstRectangularBoxGeometry & box_,
// 			   const InterpolationInfo::Method & method_ ,
// 			   const InterpolationInfo::Order  & order_ ,
// 			   std::vector<PositionParticle * >::const_iterator begin_,
// 			   std::vector<PositionParticle * >::const_iterator end_)
// {
//   begin = begin_; end = end_; box=&box_; method=method_; rcut=rcut_; 
//   K=K_; order=order_;
  
// //   std::vector<std::vector<double > > baseVector;
// //   box_.getBaseVector (baseVector);
// //   box->getBaseVector (baseVector);  
  
//   std::vector<StandardParticle * > partPool;
//   for (std::vector<PositionParticle * >::const_iterator p = begin_; 
//        p != end_; ++ p){
//     partPool.push_back (static_cast<StandardParticle * > (*p));
//   }
//   Compare com;
//   std::sort (partPool.begin(), partPool.end(), com);
//   std::vector<StandardParticle * >::iterator newend = std::unique (partPool.begin(), partPool.end());
//   partPool.erase (newend, partPool.end());

//   std::vector<double > qs;
//   for (std::vector<StandardParticle * >::iterator pp = partPool.begin();
//        pp != partPool.end(); ++pp){
//     qs.push_back ((*pp)->charge());
//   }
//   std::sort (qs.begin(), qs.end());
//   c2 = 0;
//   for (std::vector<double >::iterator i = qs.begin();
//        i != qs.end(); ++ i){
//     c2 += *i * *i;
//   }
//   N = partPool.size();
// }

// inline value_type BestBetaFinder::betaError (const value_type & beta)
// {
//   ElectrostaticInteraction ele;

// //   std::vector<std::vector<double > > baseVector;
// //   box->getBaseVector (baseVector);

//   ele.init (beta, K, rcut, *box, method, order);
// //   ele.registerParticle (begin, end);
//   double dir, rec, interpolrec;
//   ele.errorEstimatePME (c2, N, dir, rec);
//   ele.errorEstimateInterpolRec (c2, N, interpolrec);
//   return dir + interpolrec;
// }

// value_type BestBetaFinder::calBeta (const value_type & initBeta)
// {
//   value_type newbeta = initBeta;
//   value_type oldbeta ;
//   value_type h = 1e-5;
//  //  value_type epsilon = 5;
  
// //   do {
// //     oldbeta = newbeta;
// //     value_type grad = (betaError(oldbeta+h) - betaError(oldbeta-h)) / h * 0.5;
// //     newbeta = oldbeta - epsilon * grad;
// //     std::cout << "newbeta is " << newbeta << '\t' 
// // 	      << fabs (newbeta - oldbeta) << std::endl;
// //   } while (fabs (newbeta - oldbeta) > 1e-3);
 
//   // Newton iteration
//   do {
//     oldbeta = newbeta;
//     value_type v2 = betaError(oldbeta+h);
//     value_type v1 = betaError(oldbeta);
//     value_type v0 = betaError(oldbeta-h);
//     value_type grad = (v2 - v0) / h * 0.5;
//     value_type hess = (v2 - 2 * v1 + v0) / h / h;
//     newbeta = oldbeta - grad / hess;
//     if (newbeta < 0){
//       newbeta = 0.5 ;
//     }
//     if (newbeta > 3){
//       newbeta = 3;
//     }
//     std::cout << "newbeta is " << newbeta << '\t' 
// 	      << fabs (newbeta - oldbeta) << std::endl;
//   } while (fabs (newbeta - oldbeta) > 1e-5);
//   return newbeta;
// }


void ElectrostaticInteraction::get_time (double & dir_time_, double & rec_time_, 
					 double & rec_loop_time_, double & rec_conv_time_) const 
{
  dir_time_ = dir_time; 
  rec_time_ = rec_time;
  rec->get_time (rec_loop_time_, rec_conv_time_);
}

double ElectrostaticInteraction::watch_resolution () 
{ 
  return watch.resolution() < rec->watch_resolution() ? 
      rec->watch_resolution() : 
      watch.resolution();
}
