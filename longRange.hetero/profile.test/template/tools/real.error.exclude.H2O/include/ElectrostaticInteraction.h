#ifndef __wanghan_ElectrostaticInteraction_h__
#define __wanghan_ElectrostaticInteraction_h__

#include "Particle.h"
#include "BoxGeometry.h"
#include <fftw3.h>
#include "Polynominal.h"
#include "Stopwatch.h"
#include "ShortRangePairFinder.h"
#include "ToolBox.h"

class ElectrostaticInteraction_rec 
{
public:
  virtual ~ElectrostaticInteraction_rec () {};
public:
  virtual value_type calPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) = 0;
  virtual void applyInteraction (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) = 0;
  virtual value_type applyInteractionCalPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) = 0;
public:
  virtual void reset_time  () {}
  virtual void reset_watch () {}
  virtual void get_time (double & loop_time_, double & conv_time_) const {}
  virtual double watch_resolution () {return 0;}
}
  ;


class ElectrostaticInteraction_rec_ES : public ElectrostaticInteraction_rec
{
  std::vector<unsigned > K;
  value_type beta;
  value_type V;
  std::vector<std::vector<value_type > > vecA;
  std::vector<std::vector<value_type > > vecAStar;
private:
  void calV();
  void calAStar();
  
public:
  void init (const std::vector<std::vector<value_type > > &vecA_, 
	     const std::vector<unsigned > K_,
	     const value_type & beta_);
  void reinit (const std::vector<std::vector<value_type > > &vecA_, 
	       const std::vector<unsigned > K_,
	       const value_type & beta_);
public:
  value_type calPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) ;
  void applyInteraction (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) ;
  value_type applyInteractionCalPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool);
}
  ;



class EleEcalculator
{
  std::vector<double > e1K0, e1K1, e1K2;
  std::vector<double > e2K0, e2K1, e2K2;
  std::vector<double > e2pK0, e2pK1, e2pK2;
  std::vector<double > e2ppK0, e2ppK1, e2ppK2;
  std::vector<int > K;
  std::vector<int > sw; // the switch of position
  int n;
public :
  void reinit (const std::vector<unsigned > & K_,
	       const unsigned & n_);
  void get_e1 (const int & m0, const int & m1, const int & m2,
	       double & e1K0_, double & e1K1_, double & e1K2_) 
      { e1K0_ = e1K0[m0 + sw[0]]; e1K1_ = e1K1[m1 + sw[1]]; e1K2_ = e1K2[m2 + sw[2]]; }
  void get_e2 (const int & m0, const int & m1, const int & m2,
	       double & e2K0_, double & e2K1_, double & e2K2_)
      { e2K0_ = e2K0[m0 + sw[0]]; e2K1_ = e2K1[m1 + sw[1]]; e2K2_ = e2K2[m2 + sw[2]]; }
  void get_e2p (const int & m0, const int & m1, const int & m2,
		double & e2pK0_, double & e2pK1_, double & e2pK2_)
      { e2pK0_ = e2pK0[m0 + sw[0]]; e2pK1_ = e2pK1[m1 + sw[1]]; e2pK2_ = e2pK2[m2 + sw[2]]; }
  void get_e2pp (const int & m0, const int & m1, const int & m2,
		 double & e2ppK0_, double & e2ppK1_, double & e2ppK2_)
      { e2ppK0_ = e2ppK0[m0 + sw[0]]; e2ppK1_ = e2ppK1[m1 + sw[1]]; e2ppK2_ = e2ppK2[m2 + sw[2]]; }
}
    ;


class ElectrostaticInteraction_rec_BSpline : public ElectrostaticInteraction_rec
{
  std::vector<unsigned > K;
  value_type beta;
  value_type V;
  std::vector<std::vector<value_type > > vecA;
  std::vector<std::vector<value_type > > vecAStar;
  smoothInterpolationBase * Mn;
  std::vector<value_type > b0;
  std::vector<value_type > b1;
  std::vector<value_type > b2;
private:
  value_type * Q;	// do not forget to set the image part to 0
  fftw_complex * QF;
  fftw_complex * psiF;
  fftw_complex * phiF0;
  fftw_complex * phiF1;
  fftw_complex * phiF2;

  fftw_complex * phiFr0;
  fftw_complex * phiFr1;
  fftw_complex * phiFr2;
  
  fftw_complex * QFProdPsiF;
  fftw_complex * QFProdPhiF0;
  fftw_complex * QFProdPhiF1;
  fftw_complex * QFProdPhiF2;
  
  value_type * QConvPsi;
  value_type * QConvPhi0;
  value_type * QConvPhi1;
  value_type * QConvPhi2;
  
  fftw_plan forwardQ;

  fftw_plan backwardQFProdPsiF;
  fftw_plan backwardQFProdPhiF0;
  fftw_plan backwardQFProdPhiF1;
  fftw_plan backwardQFProdPhiF2;

  fftw_plan backward_phiF0;
  fftw_plan backward_phiF1;
  fftw_plan backward_phiF2;

  bool unBuild;

  Stopwatch watch;
  double conv_time;
  double loop_time;

  EleEcalculator ecal;
private:
  // values prepared during initialization
  void calV();
  void calAStar();
  void calPsiFPhiF();
  // calculated at each step. the image parts that are setted to 0 are
  // not setted again in this function
  void calQ(
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool);
  void calB();
private:
  double fenzi (const value_type & u, const unsigned & _2m);
  double fenmu (const value_type & u, const unsigned & _2m);
  double avg_e (const value_type & u, const unsigned & _2m);
  double avg_e2 (const value_type & u, const unsigned & _2m);
///////////////////////////////////////////////////////////////////
public:
  ElectrostaticInteraction_rec_BSpline () {};
  ~ElectrostaticInteraction_rec_BSpline () { clear(); }
public:
  void init (const std::vector<std::vector<value_type > > &vecA_, 
	     const std::vector<unsigned > K_,
	     const value_type & beta_,
	     const InterpolationInfo::Order & interpOrder_);
  void clear ();
  void reinit (const std::vector<std::vector<value_type > > &vecA_, 
	       const std::vector<unsigned > K_,
	       const value_type & beta_,
	       const InterpolationInfo::Order & interpOrder_) 
      {clear(); init(vecA_, K_, beta_, interpOrder_);}
public:
  void build ();
public:
  value_type calPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) ;
  void applyInteraction (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) ;
  value_type applyInteractionCalPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool);
public:
  value_type errorEstimate (
      std::vector<StandardParticle * >::const_iterator beginPool,
      std::vector<StandardParticle * >::const_iterator endPool);
  value_type errorEstimate (const double & c2, const int & N);
//   value_type errorEstimateDir (double rcut);
public:
  void reset_time  () {conv_time = 0; loop_time = 0; }
  void reset_watch () {conv_time = 0; loop_time = 0; }
  void get_time (double & loop_time_, double & conv_time_) const 
      { conv_time_ = conv_time; loop_time_ = loop_time;}
  double watch_resolution () {return watch.resolution();}
}
    ;


class EpsilonP 
{
  std::vector<std::vector<double > > ep0, ep1, ep2;
  std::vector<int > K;
  std::vector<int > sw;
  std::vector<double > h ;
  std::vector<double > hi;
  int n;
public:
  void reinit (const std::vector<unsigned > K, 
	       const unsigned & n,
	       const std::vector<unsigned > ngrind);
  void get_ep (const double & u0, const double & u1, const double & u2,
	       const int & m0, const int & m1, const int & m2,
	       double & ep0_, double & ep1_, double & ep2_) 
      { ep0_ = ep0[m0+sw[0]][int(u0*hi[0])];
	ep1_ = ep1[m1+sw[1]][int(u1*hi[1])];
	ep2_ = ep2[m2+sw[2]][int(u2*hi[2])];}
}
    ;



class ElectrostaticInteraction_rec_FBSpline : public ElectrostaticInteraction_rec
{
  std::vector<unsigned > K;
  value_type beta;
  value_type V;
  std::vector<std::vector<value_type > > vecA;
  std::vector<std::vector<value_type > > vecAStar;
  smoothInterpolationBase * Mn;
  // BSpline8 Mn;
  std::vector<value_type > b0;
  std::vector<value_type > b1;
  std::vector<value_type > b2;
private:
  value_type * Q;	// do not forget to set the image part to 0
  fftw_complex * QF;
  fftw_complex * psiF;
  fftw_complex * QFProdPsiF;
  value_type * QConvPsi;
  
  fftw_plan forwardQ;
  fftw_plan backwardQFProdPsiF;

  bool unBuild;

  Stopwatch watch;
  double conv_time;
  double loop_time;

  EleEcalculator ecal;
  EpsilonP epcal;
private:
  // values prepared during initialization
  void calV();
  void calAStar();
  void calPsiFPhiF();
  // calculated at each step. the image parts that are setted to 0 are
  // not setted again in this function
  void calQ(
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool);
  void calB();

  double gradE2 (std::vector<double > m, std::vector<double > K, double n, std::vector<double > r);
  void gradEgrade (std::vector<double > m, std::vector<double > K, double n, std::vector<double > r, double & s, double & b, double & lala, double &ahha);
  void selfTerm(std::vector<double > m, std::vector<double > K, double n, std::vector<double > r, std::vector<double > & resultre, std::vector<double > & resultim);

private:
  std::vector<double > error (const value_type & u, const value_type & x, const unsigned & _2m);
  std::vector<double > derror (const value_type & u, const value_type & x, const unsigned & _2m);
  double fenzi (const value_type & u, const unsigned & _2m);
  double fenmu (const value_type & u, const unsigned & _2m);
  double avg_e (const value_type & u, const unsigned & _2m);
  double avg_e2 (const value_type & u, const unsigned & _2m);
  double avg_epe (const value_type & u, const unsigned & _2m);
  double avg_epep (const value_type & u, const unsigned & _2m);
//////////////////////////////////////////////////////////////////////////
public:
  ElectrostaticInteraction_rec_FBSpline () {};
  ~ElectrostaticInteraction_rec_FBSpline () {clear();}
public:
  void init (const std::vector<std::vector<value_type > > &vecA_, 
	     const std::vector<unsigned > K_,
	     const value_type & beta_,
	     const InterpolationInfo::Order & interpOrder_);
  void clear ();
  void reinit (const std::vector<std::vector<value_type > > &vecA_, 
	       const std::vector<unsigned > K_,
	       const value_type & beta_,
	       const InterpolationInfo::Order & interpOrder_) 
      {clear(); init(vecA_, K_, beta_, interpOrder_);}
public:
  void build ();
public:
  value_type calPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) ;
  void applyInteraction (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool) ;
  value_type applyInteractionCalPotential (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool);
public:
  value_type calPotential_tmp (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool);
  void calForce_tmp (
      std::vector<StandardParticle * >::iterator beginPool,
      std::vector<StandardParticle * >::iterator endPool,
      std::vector<double > & force);
public:
  value_type errorEstimate (
      std::vector<StandardParticle * >::const_iterator beginPool,
      std::vector<StandardParticle * >::const_iterator endPool);
  value_type errorEstimate (const value_type & c2, const int & N);
//   value_type errorEstimateDir (double rcut);
public:
  void reset_time  () {conv_time = 0; loop_time = 0; }
  void reset_watch () {conv_time = 0; loop_time = 0; }
  void get_time (double & loop_time_, double & conv_time_) const 
      { conv_time_ = conv_time; loop_time_ = loop_time;}
  double watch_resolution () {return watch.resolution();}
}
    ;


class ElectrostaticInteraction_dir
{
  value_type rcut;
  value_type beta;
  value_type rcut2;
  value_type sqrtPIi;
public:
  ElectrostaticInteraction_dir (const value_type & rcut_,
				const value_type & beta_);
  void applyInteraction (StandardParticle & p0, 
			 StandardParticle & p1, 
			 const TimeConstRectangularBoxGeometry & boxg);
  value_type calPotential (const StandardParticle & p0, 
			   const StandardParticle & p1,
			   const TimeConstRectangularBoxGeometry & boxg);
  value_type applyInteractionCalPotential (StandardParticle & p0, 
					   StandardParticle & p1, 
					   const TimeConstRectangularBoxGeometry & boxg);
};


class ElectrostaticInteraction 
{
  value_type beta;
  std::vector<unsigned > K;
  value_type rcut;
  const TimeConstRectangularBoxGeometry * box;
  std::vector<std::vector<value_type > > baseVector;
  std::vector<std::vector<value_type > > vecAStar;
  double V;
  InterpolationInfo::Method method;
  InterpolationInfo::Order order;
  void calV();
  void calAStar();
private:
  Stopwatch watch;
  double rec_time;
  double dir_time;
  double corr_time;
  class integralTarget {
public:
    value_type beta;
    std::vector<std::vector<value_type > > vecAStar;
    integralTarget (ElectrostaticInteraction & ele);
    double operator()(double m0, double m1, double m2) const;
  };
  friend class integralTarget;
private:
  ElectrostaticInteraction_dir * dir;
  ElectrostaticInteraction_rec * rec;
  NeighborListPairFinder * pairFinder;
  std::vector<StandardParticle * > partPool;
public:
  ElectrostaticInteraction () ;
  ~ElectrostaticInteraction () {clear();}
public:
  void init (const value_type & beta_,
	     const std::vector<unsigned > & K,
	     const value_type & rcut, 
	     const TimeConstRectangularBoxGeometry & box,
	     const InterpolationInfo::Method & method = InterpolationInfo::BSpline,
	     const InterpolationInfo::Order  & order = InterpolationInfo::Sixth);
  void clear () ;
  void reinit (const value_type & beta_,
	       const std::vector<unsigned > & K,
	       const value_type & rcut, 
	       const TimeConstRectangularBoxGeometry & box,
	       const InterpolationInfo::Method & method = InterpolationInfo::BSpline,
	       const InterpolationInfo::Order  & order = InterpolationInfo::Sixth)
      {clear(); init (beta_, K, rcut, box, method, order);}
  void registerParticle (StandardParticle & p) ;
  void unregisterParticle (StandardParticle & p) ;
public:
  void build ();
  void buildRec ();
  void buildDir ();
  void rebuildDir ();
public:
  value_type calPotential ();
  void applyInteraction ();
  value_type applyInteractionCalPotential ();
  value_type applyInteractionCalPotential (std::vector<std::vector<double > > & dir_force,
					   std::vector<std::vector<double > > & rec_force);
public:
  void errorEstimatePME (const value_type & c2, const value_type & N, 
			 value_type & errorDir, value_type & errorRec);
  void errorEstimatePME (value_type & errorDir,
			 value_type & errorRec);
  void errorEstimateInterpolRec (const value_type & c2, const value_type & N,
				 value_type & errorInterpolRec);
  void errorEstimateInterpolRec (value_type & errorInterpolRec);
public:
  void reset_time  () {rec->reset_time(); rec_time = 0; dir_time = 0; corr_time = 0;}
  void reset_watch () {this->reset_time();}
  void get_time (double & dir_time_, double & rec_time_, 
		 double & rec_loop_time_, double & rec_conv_time_) const ;
  double watch_resolution () ;
};



class TimeCounter
{
  ElectrostaticInteraction ele;
  const TimeConstRectangularBoxGeometry * box;
  InterpolationInfo::Method method;
  std::vector<StandardParticle * > partPool;
private:
  value_type rcut;
  std::vector<unsigned > K;
  InterpolationInfo::Order order;
  value_type beta;
private:
  double dirtime;
  double looptime;
  double convtime;
private:
  double lower;
  double upper;
  int Klimit;
  int times;
private:
  void calTime (const value_type & myrcut, 
		const std::vector<unsigned > & myK,
		const InterpolationInfo::Order  & myorder,
		const int & mytimes,
		double & dirtime ,
		double & looptime, 
		double & convtime);
  void printStatus (double dirtime, double looptime, double convtime);
public :
  // box and method are treated as constant during calculation.
  void init (const std::vector<unsigned > & K,
	     const value_type & rcut, 
	     const TimeConstRectangularBoxGeometry & box,
	     const InterpolationInfo::Method & method = InterpolationInfo::BSpline,
	     const InterpolationInfo::Order  & order = InterpolationInfo::Sixth);
  void registerParticle (StandardParticle & p);
  void build ();
public:
  double cal_time (const value_type & myrcut, 
		   const std::vector<unsigned > & myK,
		   const InterpolationInfo::Order  & myorder ,
		   double & mydirtime,
		   double & mylooptime,
		   double & myconvtime);
}
    ;


class EleOptimizer
{
  ElectrostaticInteraction ele;
  TimeCounter tc;
  ToolBox::IntegerDecomposer id;
private:
  const TimeConstRectangularBoxGeometry * box;
  InterpolationInfo::Method method;
  std::vector<StandardParticle * > partPool;
private:
  value_type moderateRcut;
  std::vector<unsigned > moderateK;
  InterpolationInfo::Order moderateOrder;
  value_type moderateBeta;
private:
  value_type c2;
  int N;
  value_type maxrcut;
  value_type minbeta;
  value_type maxbeta;
  value_type beta2Rcut_globalBeta;
  value_type maxKload;
  value_type higherRequire;
  value_type lowerRequire;
  value_type tolScale;
  value_type timeDiffScaleLower ;
public:
  void findMinBeta ();
//   value_type beta2RcutTarget (value_type  myrcut);
//   int beta2RcutJudge (value_type  error);
  value_type beta2Rcut (const value_type & beta);
  struct beta2RcutJudge
  {
    EleOptimizer & myoptimizer;
    beta2RcutJudge (EleOptimizer & o) : myoptimizer(o) {};
    int operator () (value_type error);
  };
  friend class beta2RcutTarget;
public:
  bool beta2KO (const value_type & beta, 
		std::vector<unsigned > & K,
		InterpolationInfo::Order & order);
  bool load2K (const value_type & load,
	       const std::vector<value_type > & boxlratio,
	       std::vector<unsigned > & K);
  void bisectK (std::vector<unsigned > lowerK,
		std::vector<unsigned > upperK,
		const InterpolationInfo::Order & order,
		const value_type & beta,
		const std::vector<value_type > & boxlratio,
		std::vector<unsigned > & resultK );
  struct ParamUnit;
  struct optiParamJudge 
  {
    EleOptimizer & myoptimizer;
    double dirtime, looptime, convtime;
    value_type  rcut;
    std::vector<unsigned >  K;
    InterpolationInfo::Order  order;
    optiParamJudge (EleOptimizer & o) : myoptimizer(o) {};
    int judge (value_type beta);
    double operator () (value_type beta);
    double operator () (ParamUnit & beta);
  };
  struct ParamUnit
  {
    value_type rcut;
    std::vector<unsigned > K;
    InterpolationInfo::Order order;
    value_type beta;
    void copyFormJudge (optiParamJudge & judge) {
      rcut = judge.rcut;
      K = judge.K;
      order = judge.order;
    }
  };
  void method0618 (ParamUnit min, double tmin,  
		   ParamUnit max, double tmax,
		   optiParamJudge & judge,
		   ParamUnit & result, 
		   value_type tol);
public :
  // box and method are treated as constant during calculation.
  void init (const TimeConstRectangularBoxGeometry & box_,
	     const InterpolationInfo::Method & method_ ,
	     const value_type & requirePrecision,
	     const std::vector<unsigned > & moderateK_ = std::vector<unsigned > (3, 64),
	     const InterpolationInfo::Order & moderateOrder_ = 6,
	     const value_type & moderateRcut_ = 3,
	     const std::vector<unsigned >  & maxK = std::vector<unsigned > (3, 300),
	     const InterpolationInfo::Order & maxOrder_ = 8,
	     const value_type & maxRcut_ = 10,
	     const value_type & maxBeta_ = 2);
  void registerParticle (StandardParticle & p);
  void build ();
  bool optiParam (value_type & rcut,
		  std::vector<unsigned > & K,
		  InterpolationInfo::Order & order, 
		  value_type & beta,
		  const value_type & betaTol = 0.01);
}
    ;


  
      

// class BestBetaFinder 
// {
//   double c2;
//   int N;
//   std::vector<PositionParticle * >::const_iterator begin;
//   std::vector<PositionParticle * >::const_iterator end;
//   const TimeConstRectangularBoxGeometry * box;
//   InterpolationInfo::Method method;
//   value_type rcut;
//   std::vector<unsigned > K;
//   InterpolationInfo::Order order;
// public:
//   void init (const std::vector<unsigned > & K_,
// 	     const value_type & rcut_, 
// 	     const TimeConstRectangularBoxGeometry & box_,
// 	     const InterpolationInfo::Method & method_,
// 	     const InterpolationInfo::Order  & order_,
// 	     std::vector<PositionParticle * >::const_iterator begin_,
// 	     std::vector<PositionParticle * >::const_iterator end_);
//   inline value_type betaError (const value_type & beta);
//   value_type calBeta (const value_type & initBeta);
// }
//     ;


// class KOFinder // find mesh size and order of interpolation
// {
// public:
//   void init (const value_type & rcut_, 
// 	     const TriclinicBoxGeometry & box_,
// 	     const InterpolationInfo::Method & method_,
// 	     std::vector<PositionParticle * >::const_iterator begin_,
// 	     std::vector<PositionParticle * >::const_iterator end_);
// }
	     
  



  
#endif 
  



// class ElectrostaticInteraction_rec_Lagrange : public ElectrostaticInteraction_rec 
// {
//   std::vector<unsigned > K;
//   value_type beta;
//   value_type V;
//   std::vector<std::vector<value_type > > vecA;
//   std::vector<std::vector<value_type > > vecAStar;
// private:
//   Polynominalp3 W2p;
// private:
//   value_type * Q;	// do not forget to set the image part to 0
//   fftw_complex * QF;
//   fftw_complex * psiF;
//   fftw_complex * phiF0;
//   fftw_complex * phiF1;
//   fftw_complex * phiF2;
  
//   fftw_complex * QFProdPsiF;
//   fftw_complex * QFProdPhiF0;
//   fftw_complex * QFProdPhiF1;
//   fftw_complex * QFProdPhiF2;
  
//   value_type * QConvPsi;
//   value_type * QConvPhi0;
//   value_type * QConvPhi1;
//   value_type * QConvPhi2;
  
//   fftw_plan forwardQ;

//   fftw_plan backwardQFProdPsiF;
//   fftw_plan backwardQFProdPhiF0;
//   fftw_plan backwardQFProdPhiF1;
//   fftw_plan backwardQFProdPhiF2;
// private:
//   // values prepared during initialization
//   void calV();
//   void calAStar();
//   void calPsiFPhiF();
//   // calculated at each step. the image parts that are setted to 0 are
//   // not setted again in this function
//   void calQ();
// public:
//   void init (const std::vector<std::vector<value_type > > &vecA_, 
// 	     const std::vector<unsigned > K_,
// 	     const value_type & beta_);
//   void clear ();
//   void reinit (const std::vector<std::vector<value_type > > &vecA_, 
// 	       const std::vector<unsigned > K_,
// 	       const value_type & beta_) 
//       {clear(); init(vecA_, K_, beta_);}
// public:
//   ElectrostaticInteraction_rec_Lagrange () {};
//   ~ElectrostaticInteraction_rec_Lagrange () {clear();}
// public :
//   void applyInteraction (const double & time, const TriclinicBoxGeometry & box);
//   value_type calPotential (const double & time, const TriclinicBoxGeometry & box);
//   value_type calPotential_tmp (const double & time, const TriclinicBoxGeometry & box);
//   void calForce_tmp (const double & time, const TriclinicBoxGeometry & box, std::vector<double > & force);
// public:
//   value_type errorEstimate ();
//   value_type errorEstimateDir (double rcut);
// }
//     ;


