#include "ElectrostaticInteraction.h"
#include "ToolBox.h"
#include <cmath>

void TimeCounter::init (
    const std::vector<unsigned > & K_,
    const value_type & rcut_, 
    const TimeConstRectangularBoxGeometry & box_,
    const InterpolationInfo::Method & method_,
    const InterpolationInfo::Order  & order_ )
{
  beta = 2;
  K = K_;
  rcut = rcut_;
  box = &box_;
  method = method_;
  order = order_;
  
  ele.init (beta, K, rcut, *box, method, order);
  
  lower = ele.watch_resolution() * 50;
  upper = lower * 30;
  times = 1;
  Klimit = 90;
}

void TimeCounter::registerParticle (StandardParticle & p)
{
  partPool.push_back (&p);
}

void TimeCounter::calTime (
    const value_type & myrcut, 
    const std::vector<unsigned > & myK,
    const InterpolationInfo::Order  & myorder,
    const int & mytimes,
    double & dirtime ,
    double & looptime, 
    double & convtime)
{
  ele.reinit (beta, myK, myrcut, *box, method, myorder);
  for (std::vector<StandardParticle * >::iterator i = partPool.begin();
       i != partPool.end(); ++i){
    ele.registerParticle (**i);
  }
  ele.build();
  
  ele.reset_watch();
  for (int i = 0; i < mytimes; ++i){
    ele.applyInteraction();
  }
  double rectime;
  ele.get_time (dirtime, rectime, looptime, convtime);
}

void TimeCounter::printStatus (double dirtime, double looptime, double convtime)
{  
#ifdef VERBOSE_INFO
  std::cout.precision (3);
  std::cout << "times " << times
	    << "\trcut " << rcut 
	    << "\torder " << order 
	    << "\tK " << K[0] << std::endl;
  std::cout << "now, the statue is \t" ;
  std::cout << dirtime << "\t"
	    << looptime << "\t"
	    << convtime << std::endl;
  std::cout.precision (16);
#endif
}


void TimeCounter::build ()
{
  ToolBox::IntegerDecomposer id;
  std::vector<int > index;
  std::vector<value_type > boxl;
  box->getBoxSize (boxl);
  std::sort (boxl.begin(), boxl.end());
  bool rcutUpperNo = false;
  bool KUpperNo = false;
  bool moreOpt = true;
#ifdef  VERBOSE_INFO
  std::cout << "upper time is " << upper
	    << "\nlower time is " << lower << std::endl;
#endif  
  calTime (rcut, K, order,  1, dirtime, looptime, convtime);
  printStatus (dirtime, looptime, convtime);
  
  
  while (moreOpt){
    moreOpt = false;
    if (dirtime > upper){
      rcut /= 2;
      moreOpt = true;
    }
    if (looptime > upper) {
      if (order > 2){
	order -= 2;
      }
      moreOpt = true;
    }
    if (convtime > upper){
      for (int i = 0; i < 3; ++i){
	if (K[i] >= 16){
	  K[i] /= 2;
	}
	K[i] = id.decompose (K[i], index);
      }
      moreOpt = true;
    }
    if (dirtime < lower){
      if ((rcut * 2) < 0.9 * 0.5 * boxl[0]) {
	rcut *= 2;
	moreOpt = true;
      }
      else {
	rcut = 0.9 * 0.5 * boxl[0];
	rcutUpperNo = true;
      }
    }
    if (looptime < lower){
      if (order < 10){
	order += 2;
      }
      moreOpt = true;
    }
    if (convtime < lower){
      for (int i = 0; i < 3; ++i){
	if (K[i] * 2 > unsigned(Klimit)){
	  K[i] = Klimit;
	  K[i] = id.decompose (K[i], index);
	  KUpperNo = true;
	}
	else{
	  K[i] *= 2;
	  K[i] = id.decompose (K[i], index);
	  moreOpt = true;
	}
      }
    }
    if (KUpperNo || rcutUpperNo){
      times *= 2;
      KUpperNo = false;
      rcutUpperNo = false;
      moreOpt = true;
    }
    if (moreOpt){
      calTime (rcut, K, order, times, dirtime, looptime, convtime);
      printStatus (dirtime, looptime, convtime);
    }
  }

  double tmp1 = double(K[0]) * double (K[1]) * double(K[2]);
  std::cout << "## dir time constant is " 
	    << dirtime / double (times) / (rcut*rcut*rcut) << std::endl;
  std::cout << "## loop time constant is " 
	    << looptime / double (times) / double(order) / double(order) / double(order)  << std::endl;
  std::cout << "## conv time constant is "
	    << convtime / double (times) / tmp1 / log(tmp1) << std::endl;
}

double TimeCounter::cal_time (const value_type & myrcut, 
			      const std::vector<unsigned > & myK,
			      const InterpolationInfo::Order  & myorder ,
			      double & mydirtime,
			      double & mylooptime,
			      double & myconvtime)
{
//   std::cout << looptime << std::endl;
  mydirtime = 
      (myrcut / rcut) * 
      (myrcut / rcut) * 
      (myrcut / rcut) * 
      dirtime / double (times);
  mylooptime = 
      (double(myorder) / double(order)) * 
      (double(myorder) / double(order)) * 
      (double(myorder) / double(order)) * 
      looptime / double (times);
  double tmp1 = double(K[0]) * double (K[1]) * double(K[2]);
  double tmp2 = double(myK[0]) * double (myK[1]) * double(myK[2]);
  myconvtime = tmp2 / tmp1 * log(tmp2) / log(tmp1) *
      convtime / double (times);
  return mydirtime + mylooptime + myconvtime;
}

