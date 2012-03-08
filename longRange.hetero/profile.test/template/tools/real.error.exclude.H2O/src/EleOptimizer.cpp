#include "ElectrostaticInteraction.h"
#include "VectorOperation.h"

void EleOptimizer::init (const TimeConstRectangularBoxGeometry & box_,
			 const InterpolationInfo::Method & method_ ,
			 const value_type & requirePrecision,
			 const std::vector<unsigned > & moderateK_,
			 const InterpolationInfo::Order & moderateOrder_,
			 const value_type & moderateRcut_,
			 const std::vector<unsigned >  & maxK_,
			 const InterpolationInfo::Order & maxOrder_,
			 const value_type & maxRcut_,
			 const value_type & maxBeta_)
{
  moderateK = moderateK_;
  moderateRcut = moderateRcut_;
  box = &box_;
  method = method_;
  moderateOrder = moderateOrder_;
//   std::cout << maxBeta_ << "::::::::::::::::::::::::::::::::::::::" << std::endl;
  maxbeta = maxBeta_;
  maxKload = double(maxK_[0]) * double (maxK_[1]) * double(maxK_[2]);
  
  timeDiffScaleLower = 0.9;
  
  std::vector<value_type > boxl;
  box->getBoxSize (boxl);
  std::sort (boxl.begin(), boxl.end());
  maxrcut = (0.45 * boxl[0] < maxRcut_) ?  (0.45 * boxl[0]) : (maxRcut_);
  
  ele.init (2, moderateK, moderateRcut, *box, method, moderateOrder);
  tc.init  (   moderateK, moderateRcut, *box, method, moderateOrder);
  
  higherRequire = requirePrecision * 0.5;
  tolScale = 0.1;
  lowerRequire = higherRequire - tolScale * higherRequire;
}

void EleOptimizer::registerParticle (StandardParticle & p)
{
  partPool.push_back (&p);
  tc.registerParticle (p);
}

void EleOptimizer::build ()
{
#ifdef VERBOSE_INFO
  std::cout << "### TC is building" << std::endl;
#endif
  tc.build();
#ifdef VERBOSE_INFO
  std::cout << "### finished " << std::endl;
#endif
  c2 = 0;
  for (std::vector<StandardParticle *>::iterator ppart = partPool.begin();
       ppart != partPool.end(); ++ppart){
    c2 += (*ppart)->charge() * (*ppart)->charge() ;
  }
  N = partPool.size();
  
  findMinBeta ();
}


// higherRequire should be set in advance
void EleOptimizer::findMinBeta()
{
  value_type mybeta = 0.4;
  value_type direrror = 2 * higherRequire, recerror;
  
  while (direrror > higherRequire){
    mybeta += 0.05;
    if (mybeta > 4){
      std::cerr << "cannot find a minimal beta, exit" << std::endl;
      exit (1); 
    }
    ele.reinit (mybeta, moderateK, maxrcut, *box, method, moderateOrder);
    ele.errorEstimatePME (c2, N, direrror, recerror);
  }
#ifdef VERBOSE_INFO
  std::cout << "at minimal beta " << mybeta
	    << " the dire error is " << direrror << std::endl;
  #endif
  moderateBeta = 0.5 * (mybeta + maxbeta);
  minbeta = mybeta;
}


int EleOptimizer::beta2RcutJudge::operator () (value_type  myrcut)
{
  myoptimizer.ele.reinit (myoptimizer.beta2Rcut_globalBeta, 
			  myoptimizer.moderateK, 
			  myrcut, 
			  *(myoptimizer.box), 
			  myoptimizer.method, 
			  myoptimizer.moderateOrder);
  double error, b;
  myoptimizer.ele.errorEstimatePME (myoptimizer.c2, myoptimizer.N, error, b);

  if ((error - myoptimizer.lowerRequire) > 
      myoptimizer.tolScale * myoptimizer.higherRequire)
    return 1;
  else if ((error - myoptimizer.lowerRequire) < 
	   - myoptimizer.tolScale * myoptimizer.higherRequire)
    return -1;
  else 
    return 0;
}
  

template<typename D2I>
value_type bisect (value_type min, value_type max,
		   D2I judge)
{
  int c;
  int count = 0;
  while (1){
    value_type mid = 0.5 * (min + max);
    if ((c=judge(mid)) == 0) 
      return mid;
    else if (c == -1) 
      max = mid;
    else 
      min = mid;
    if (++count > 50){
      std::cerr << "cannot find a solution" << std::endl;
      return mid;
    }
  }
}

// template<typename D2I>
// value_type bisect1 (value_type min, value_type max,
// 		    D2I judge)
// {
//   int c;
//   int count = 0;
//   while (1){
//     value_type mid = 0.5 * (min + max);
//     c=judge(mid);
//     std::cout << "rcut " << judge.rcut
// 	      << "\torder " << judge.order
// 	      << "\tK " << judge.K[0]  << " " << judge.K[1] << " " << judge.K[2]
// 	      << "\tdir " << judge.dirtime 
// 	      << "\tloop " << judge.looptime 
// 	      << "\tconv " << judge.convtime 
// 	      << "\trec " << judge.looptime + judge.convtime 
// 	      << "\ttotal " << judge.dirtime + judge.looptime + judge.convtime 
// 	      << std::endl;
//     if (c == 0) 
//       return mid;
//     else if (c == -1) 
//       max = mid;
//     else 
//       min = mid;
//     if (++count > 50){
//       std::cerr << "cannot find a solution" << std::endl;
//       return mid;
//     }
//   }
// }


value_type EleOptimizer::beta2Rcut (const value_type & beta)
{
//   findMinBeta ();
  
  beta2Rcut_globalBeta = beta;
  beta2RcutJudge  judge (*this);
  value_type minrcut = moderateRcut * 0.05;
#ifdef VERBOSE_INFO
  std::cout << "given beta is " << beta << std::endl;
  std::cout << "starting interval is " << minrcut << "\t" << maxrcut << std::endl;
  std::cout << "precision is " << lowerRequire << '\t' << tolScale << std::endl;
#endif
  if (judge(minrcut) != 1){
    return minrcut;
  }
  else if (judge(maxrcut) == 0)
    return maxrcut;
  else
    return bisect<beta2RcutJudge> (minrcut, maxrcut, judge);
}

bool EleOptimizer::load2K (const value_type & load,
			   const std::vector<value_type > & boxlratio,
			   std::vector<unsigned > & K)
{
  K.resize(3);
  value_type multiratio = boxlratio[0] * boxlratio[1] * boxlratio[2];
  value_type unitK = pow( load / multiratio, 1./3.);
  K[0] = id.decompose(int(unitK * boxlratio[0]));
  K[1] = id.decompose(int(unitK * boxlratio[1]));
  K[2] = id.decompose(int(unitK * boxlratio[2]));
  if (double(K[0]) * double(K[1]) * double(K[2]) >= maxKload) {
    unitK = pow( maxKload 
		 / multiratio
		 , 1./3.);
    K[0] = id.decompose(int(unitK * boxlratio[0]));
    K[1] = id.decompose(int(unitK * boxlratio[1]));
    K[2] = id.decompose(int(unitK * boxlratio[2]));
    return false;
  }
  return true;
}

// void EleOptimizer::load2K2 (const value_type & load,
// 			    const std::vector<value_type > & boxlratio,
// 			    std::vector<unsigned > & K)
// {
//   K.resize(3);
//   value_type multiratio = boxlratio[0] * boxlratio[1] * boxlratio[2];
//   value_type unitK = pow( load / multiratio, 1./3.);
//   K[0] = id.decompose(int(unitK * boxlratio[0]));
//   K[1] = id.decompose(int(unitK * boxlratio[1]));
//   K[2] = id.decompose(int(unitK * boxlratio[2]));
// }

void EleOptimizer::bisectK (std::vector<unsigned > lowerK,
			    std::vector<unsigned > upperK,
			    const InterpolationInfo::Order & order,
			    const value_type & beta,
			    const std::vector<value_type > & boxlratio,
			    std::vector<unsigned > & resultK )
{
#ifdef VERBOSE_INFO
  std::cout << "now in the bsect K " << std::endl;
#endif
  resultK.resize(3);
  for (int i = 0; i < 10; ++i){
    double tmp0 = lowerK[0] + upperK[0];
    double tmp1 = lowerK[1] + upperK[1];
    double tmp2 = lowerK[2] + upperK[2];
    double error;
    double dirtime, looptime, convtime;
    load2K(0.125 * tmp0 * tmp1 * tmp2, boxlratio, resultK);
    if (resultK[0] == upperK[0] || resultK[0] == lowerK[0]) {
#ifdef VERBOSE_INFO
      std::cout << "haha upper " << upperK[0] << "\tlower " << lowerK[0] << std::endl;
      std::cout << "error " << error
		<< "\torder " << order
		<< "\tK " << resultK[0]  << " " << resultK[1] << " " << resultK[2]
		<< "\tlooptime " << looptime 
		<< "\tconvtime " << convtime 
		<< "\ttotaltime " << looptime + convtime << std::endl;
#endif
      resultK = upperK; 
      return;
    }
    ele.reinit (beta, resultK, moderateRcut, 
		*box, method, order);
    ele.errorEstimateInterpolRec (c2, N, error);
    tc.cal_time (moderateRcut, resultK, order, dirtime, looptime, convtime);
#ifdef VERBOSE_INFO
    std::cout << "error " << error
	      << "\torder " << order
	      << "\tK " << resultK[0]  << " " << resultK[1] << " " << resultK[2]
	      << "\tlooptime " << looptime 
	      << "\tconvtime " << convtime 
	      << "\ttotaltime " << looptime + convtime << std::endl;
#endif
    if (error <= higherRequire && error >= higherRequire * (1-tolScale*2))
      return;
    else if (error > higherRequire){
      lowerK[0] = resultK[0];
      lowerK[1] = resultK[1];
      lowerK[2] = resultK[2];
    }
    else if (error < higherRequire * (1-tolScale*2)){
      upperK[0] = resultK[0];
      upperK[1] = resultK[1];
      upperK[2] = resultK[2];
    }
  }
}
    

bool EleOptimizer::beta2KO (const value_type & beta, 
			    std::vector<unsigned > & K,
			    InterpolationInfo::Order & order)
{
#ifdef VERBOSE_INFO
  std::cout.precision(3);
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
  std::cout << "@@@ beta to K and order, beta is " << beta << std::endl;
#endif
  std::vector<value_type > boxlratio;
  box->getBoxSize (boxlratio);
  std::vector<int > orders;
  orders.push_back(8);
  orders.push_back(6);
  orders.push_back(4);
//  orders.push_back(2);
  std::vector<unsigned > optiK;
  int optiOrder;
  double optiTime = 1e10;
  double dirtime, looptime, convtime;
  
  double tmperror;
  for (int i = 0; i < int(orders.size()); ++i){
    std::vector<unsigned > myK;
    int myorder;
    int thisorder = orders[i];
    std::vector<unsigned > lowerK ;
    std::vector<unsigned > upperK ;
    double lowerload = 4*4*4;
    double upperload = 4*4*4;
    // cal lower error
    load2K (lowerload, boxlratio, lowerK);
    upperK = lowerK;
    ele.reinit (beta, lowerK, moderateRcut, *box, method, thisorder);
    ele.errorEstimateInterpolRec (c2, N, tmperror);
    // if lower error is small enough, then use lower K and lower order
    if (tmperror <= higherRequire) {
      optiK = lowerK;
      optiOrder = thisorder;
      break;
    }
    do {
      // set lower K to be upper K
      lowerload = upperload;
      lowerK = upperK;
      // calculate upper K
      upperload = upperload * 8;
      bool hitmax;
      hitmax = ! load2K (upperload, boxlratio, upperK);
      std::cout << ";; order " << thisorder
		<< ", lower K is " << lowerK[0] 
		<< ", upper K is " << upperK[0];
      // check time of lower K, if the time is bigger than the
      // optimized one, to next order
      tc.cal_time (moderateRcut, lowerK, thisorder, dirtime, looptime, convtime);
      if (looptime + convtime > optiTime){
	std::cout << ", go to next order" << std::endl;
	goto nextorder;
      }
      // cal upper error
      ele.reinit (beta, upperK, moderateRcut, *box, method, thisorder);
      ele.errorEstimateInterpolRec (c2, N, tmperror);
      std::cout << ", upper error is " << tmperror
		<< " (required is " << higherRequire 
		<< "), hitmax " << hitmax << std::endl;
      if (i == 0 && hitmax && tmperror >= higherRequire) {
	return false;
      }
      if (i != 0 && hitmax && tmperror >= higherRequire) {
	K = optiK;
	order = optiOrder;
	return true;
      }
      // if upper error is bigger than the required error, then
      // increase upper K
    } while (tmperror >= higherRequire);
    // find solution
    bisectK (lowerK, upperK, thisorder, beta, boxlratio, myK);
    myorder = thisorder;
    tc.cal_time (moderateRcut, myK, myorder, dirtime, looptime, convtime);
    if (looptime + convtime < optiTime){
      optiTime = looptime + convtime;
      optiK = myK;
      optiOrder = myorder;
    }
    nextorder :
    std::cout << "going to next order " << std::endl;
  }
  
  K = optiK;
  order = optiOrder;
  return true;
}

      
      

// bool EleOptimizer::beta2KO (const value_type & beta, 
// 			    std::vector<unsigned > & K,
// 			    InterpolationInfo::Order & order)
// {
// #ifdef VERBOSE_INFO
//   std::cout.precision(3);
//   std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
//   std::cout << "@@@ beta to K and order, beta is " << beta << std::endl;
// #endif
//   int maxorder = 8;
//   std::vector<value_type > boxlratio;
//   box->getBoxSize (boxlratio);
//   VectorOperation::scale (boxlratio, 1./boxlratio[0]);
//   if (boxlratio[1] < 1)
//     VectorOperation::scale (boxlratio, 1./boxlratio[1]);
//   if (boxlratio[2] < 1)
//     VectorOperation::scale (boxlratio, 1./boxlratio[2]);

//   std::vector<unsigned > startK(3);
//   load2K (32*32*32, boxlratio, startK);

//   std::vector<std::vector<unsigned > > lowerKs;
//   std::vector<std::vector<unsigned > > upperKs;
//   std::vector<value_type > lowerKError;
//   std::vector<value_type > upperKError;
//   // calculate the K range for each order
//   for (int thisorder = 2; thisorder <= maxorder; thisorder+=2){
//     double dirtime, looptime, convtime;
//     tc.cal_time (moderateRcut, startK, thisorder, 
// 		 dirtime, looptime, convtime);
//     std::vector<unsigned > tmpK;
//     load2K ((1-timeDiffScaleLower) * 
// 	    (looptime/convtime * startK[0] * startK[1] *startK[2]), 
// 	    boxlratio, tmpK) ;
//     lowerKs.push_back (tmpK);
//   }
//   int thisorder;
//   for (thisorder = 2; thisorder <= maxorder-2; thisorder+=2){
//     double dirtime, looptime0, convtime0;
//     double looptime1, convtime1;
//     tc.cal_time (moderateRcut, startK, thisorder, 
// 		 dirtime, looptime0, convtime0);
//     tc.cal_time (moderateRcut, lowerKs[thisorder/2], thisorder+2, 
// 		 dirtime, looptime1, convtime1);
//     std::vector<unsigned > tmpK;
//     double requireConvtime = looptime1 + convtime1 - looptime0;
//     load2K ((requireConvtime / convtime0) * startK[0] * startK[1] *startK[2], 
// 	    boxlratio, tmpK) ;
//     upperKs.push_back (tmpK);
//   }
//   double dirtime, looptime, convtime;
//   tc.cal_time (moderateRcut, startK, thisorder, 
// 	       dirtime, looptime, convtime);
//   std::vector<unsigned > tmpK;
//   if (looptime/convtime * startK[0] * startK[1] *startK[2] > maxKload){
//     load2K (maxKload, boxlratio, tmpK);
//   }
//   else {
//     load2K ((1) * 
// 	  (looptime/convtime * startK[0] * startK[1] *startK[2]), 
// 	  boxlratio, tmpK) ;
//   }
//   upperKs.push_back (tmpK);
//   // calculate the error range for each K range
// #ifdef VERBOSE_INFO	    
//   std::cout << "required high error is " << higherRequire << std::endl;
// #endif  
//   for (int thisorder = 2; thisorder <= maxorder; thisorder+=2){
//     double errorInterMin, errorInterMax;
//     std::cout.precision(3);
//     ele.reinit (beta, lowerKs[thisorder/2-1], moderateRcut, 
// 		*box, method, thisorder);
//     ele.errorEstimateInterpolRec (c2, N, errorInterMax);
//     lowerKError.push_back(errorInterMax);
//     tc.cal_time (moderateRcut, lowerKs[thisorder/2-1], thisorder, 
// 		 dirtime, looptime, convtime);
// #ifdef VERBOSE_INFO
//     std::cout << "error " << errorInterMax
// 	      << "\torder " << thisorder 
// 	      << "\tK " << lowerKs[thisorder/2-1][0]  << " " << lowerKs[thisorder/2-1][1] << " " << lowerKs[thisorder/2-1][2]
// 	      << "\tlooptime " << looptime 
// 	      << "\tconvtime " << convtime 
// 	      << "\ttotaltime " << looptime + convtime << std::endl;
// #endif
//     std::cout.precision(3);
//     ele.reinit (beta, upperKs[thisorder/2-1], moderateRcut, 
// 		*box, method, thisorder);
//     ele.errorEstimateInterpolRec (c2, N, errorInterMin);
//     upperKError.push_back(errorInterMin);
//     tc.cal_time (moderateRcut, upperKs[thisorder/2-1], thisorder, 
// 		 dirtime, looptime, convtime);
// #ifdef VERBOSE_INFO
//     std::cout << "error " << errorInterMin
// 	      << "\torder " << thisorder 
// 	      << "\tK " << upperKs[thisorder/2-1][0]  << " " << upperKs[thisorder/2-1][1] << " " << upperKs[thisorder/2-1][2]
// 	      << "\tlooptime " << looptime 
// 	      << "\tconvtime " << convtime 
// 	      << "\ttotaltime " << looptime + convtime << std::endl;
// #endif
//     if (errorInterMax < higherRequire){
//       std::vector<unsigned > tmplowK;
//       load2K (4* 4 * 4, boxlratio, tmplowK);
//       bisectK (tmplowK, lowerKs[thisorder/2-1], thisorder,
// 	       beta, boxlratio, K);
//       order = thisorder;
//       return true;
//     }
//     if (errorInterMax > higherRequire && errorInterMin < higherRequire){
//       std::vector<unsigned > tmpK0;
//       std::vector<unsigned > tmpK1;
//       bisectK (lowerKs[thisorder/2-1], upperKs[thisorder/2-1], thisorder,
// 	       beta, boxlratio, tmpK0);
// #ifdef Ele_Optimizer_Exhaustive
//       if (thisorder > 2 && lowerKError[thisorder/2-1] > upperKError[thisorder/2-2]){
// 	std::cout << "$$$$$$$$$$$$$$$$$$" << std::endl;
// 	std::vector<unsigned > tmpupK;
// 	load2K (maxKload, boxlratio, tmpupK);
// 	value_type tmpuperror;
// 	ele.reinit (beta, tmpupK, moderateRcut, *box, method, thisorder-2);
// 	ele.errorEstimateInterpolRec (c2, N, tmpuperror);
// #ifdef VERBOSE_INFO
// 	std::cout << "$$$$$$$$$$$$$$$$$$$$ the lower error is " << tmpuperror << std::endl;
// #endif
// 	if (tmpuperror < higherRequire){
// 	  bisectK (upperKs[thisorder/2-2], tmpupK, thisorder-2, 
// 		   beta, boxlratio, tmpK1);
// 	}
//       }
// #endif
//       if (tmpK1.empty()){
// 	K = tmpK0;
// 	order= thisorder;
//       }
//       else {
// 	double dirtime0, looptime0, convtime0;
// 	double dirtime1, looptime1, convtime1;
// 	tc.cal_time (moderateRcut, tmpK0, thisorder, 
// 		     dirtime0, looptime0, convtime0);
// 	tc.cal_time (moderateRcut, tmpK1, thisorder-2, 
// 		     dirtime1, looptime1, convtime1);
// 	if (looptime1 + convtime1 < looptime0 + convtime0){
// 	  K = tmpK1;
// 	  order = thisorder-2;
// 	}
// 	else {
// 	  K = tmpK0;
// 	  order= thisorder;
// 	}
//       }  
//       return true;
//     }
//   }
//   return false;  
// }

   
int EleOptimizer::optiParamJudge::judge (value_type beta)
{
  rcut = myoptimizer.beta2Rcut (beta);
  if (! myoptimizer.beta2KO (beta, K, order))
    return -2;
  myoptimizer.tc.cal_time (rcut, K, order, dirtime, looptime, convtime);
  
  value_type mytol = 0.1;
  if (dirtime / (looptime + convtime) < 1-mytol)
    return -1;
  else if (dirtime / (looptime + convtime) > 1 + mytol)
    return 1;
  else 
    return 0;
}

double  EleOptimizer::optiParamJudge::operator() (value_type beta)
{
  rcut = myoptimizer.beta2Rcut (beta);
  myoptimizer.beta2KO (beta, K, order);
  myoptimizer.tc.cal_time (rcut, K, order, dirtime, looptime, convtime);
  return dirtime + looptime + convtime;
}

double  EleOptimizer::optiParamJudge::operator() (ParamUnit & p)
{
  p.rcut = myoptimizer.beta2Rcut (p.beta);
  myoptimizer.beta2KO (p.beta, p.K, p.order);
  myoptimizer.tc.cal_time (p.rcut, p.K, p.order, dirtime, looptime, convtime);
  return dirtime + looptime + convtime;
}


static void findMini (double a, double b, double c, double d,
		      int & mini1, int & mini2)
{
  std::vector<std::pair<double, int> > tmp;
  tmp.push_back(std::pair<double, int> (a, 1));
  tmp.push_back(std::pair<double, int> (b, 2));
  tmp.push_back(std::pair<double, int> (c, 3));
  tmp.push_back(std::pair<double, int> (d, 4));
  
  std::sort (tmp.begin(), tmp.end());
  mini2 = tmp[1].second;
  mini1 = tmp[0].second;
}

void EleOptimizer::method0618 (ParamUnit min, double tmin,  
			       ParamUnit max, double tmax,
			       optiParamJudge & judge,
			       ParamUnit & result, 
			       value_type tol)
{
  value_type tau1 = 0.5 * (sqrt(5) -1);
  value_type tau0 = 1 - tau1;
  ParamUnit s0, s1;
  double ts0, ts1;
  s0.beta = min.beta + tau0 * (max.beta - min.beta);
  s1.beta = min.beta + tau1 * (max.beta - min.beta);
  ts0 = judge (s0);
  ts1 = judge (s1);
#ifdef VERBOSE_INFO
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	    << "%%% time: " << tmin << '\t' << ts0 << "\t"
	    << ts1 << "\t" << tmax << "\n"
	    << "%%% beta: " << min.beta << "\t" << s0.beta << "\t" 
	    << s1.beta << "\t" << max.beta << "\n";
#endif  
  while (max.beta - min.beta > tol){
//     if (ts0 < ts1){
    int mini0, mini1;
    findMini (tmin, ts0, ts1, tmax, mini0, mini1);
    if ((mini0 == 1 || mini1 == 1) || 
	(mini0 != 4 && mini1 != 4 && ts0 < ts1)){
#ifdef VERBOSE_INFO
      std::cout << "%%% the two minimals are " << mini0 << " and " << mini1 
		<< ". choose interval I" << std::endl;
#endif  
      max = s1;
      tmax = ts1;
      s1 = s0;
      ts1 = ts0;
      s0.beta = min.beta + tau0 * (max.beta - min.beta);
      ts0 = judge (s0);
    }
    else {
#ifdef VERBOSE_INFO
      std::cout << "%%% the two minimals are " << mini0 << " and " << mini1 
		<< ". choose interval II" << std::endl;
#endif  
      min = s0;
      tmin = ts0;
      s0 = s1;
      ts0 = ts1;
      s1.beta = min.beta + tau1 * (max.beta - min.beta);
      ts1 = judge (s1);
    }
#ifdef VERBOSE_INFO    
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	      << "%%% time: " << tmin << '\t' << ts0 << "\t"
	      << ts1 << "\t" << tmax << "\n"
	      << "%%% beta: " << min.beta << "\t" << s0.beta << "\t" 
	      << s1.beta << "\t" << max.beta << "\n";
#endif
  }
  
  ts0 < ts1 ? result = s0 : result = s1;
}      
      
    

bool EleOptimizer::optiParam (value_type & rcut,
			      std::vector<unsigned > & K,
			      InterpolationInfo::Order & order, 
			      value_type & beta,
			      const value_type & betaTol)
{
  optiParamJudge opj (*this);
  
  ParamUnit min, max;
  double tmin, tmax;
  if (opj.judge (minbeta) == -2){
    return false;
  }
  min.beta = minbeta;
  min.rcut = opj.rcut;
  min.K = opj.K;
  min.order = opj.order;
  tmin = opj.dirtime + opj.looptime + opj.convtime ;
#ifdef VERBOSE_INFO
  std::cout << "rcut " << opj.rcut
	    << "\torder " << opj.order
	    << "\tK " << opj.K[0]  << " " << opj.K[1] << " " << opj.K[2]
	    << "\tdir " << opj.dirtime 
	    << "\tloop " << opj.looptime 
	    << "\tconv " << opj.convtime 
	    << "\trec " << opj.looptime + opj.convtime 
	    << "\ttotal " << opj.dirtime + opj.looptime + opj.convtime 
	    << std::endl;
#endif
//   std::cout << maxbeta << ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;" << std::endl;
  while (opj.judge (maxbeta) == -2){
    maxbeta = 0.8 * (maxbeta - minbeta) + minbeta;
    if (maxbeta - minbeta < betaTol) {
      beta = min.beta;
      rcut = min.rcut;
      K = min.K;
      order = min.order;
      return true;
    }
  }
  max.beta = maxbeta;
  max.rcut = opj.rcut;
  max.K = opj.K;
  max.order = opj.order;
  tmax = opj.dirtime + opj.looptime + opj.convtime ;
#ifdef VERBOSE_INFO
  std::cout << "rcut " << opj.rcut
	    << "\torder " << opj.order
	    << "\tK " << opj.K[0]  << " " << opj.K[1] << " " << opj.K[2]
	    << "\tdir " << opj.dirtime 
	    << "\tloop " << opj.looptime 
	    << "\tconv " << opj.convtime 
	    << "\trec " << opj.looptime + opj.convtime 
	    << "\ttotal " << opj.dirtime + opj.looptime + opj.convtime 
	    << std::endl;
#endif
  ParamUnit resultpara;
  method0618 (min, tmin, max, tmax, opj, resultpara, betaTol);
  
  rcut = resultpara.rcut;
  K = resultpara.K;
  order = resultpara.order;
  beta = resultpara.beta;

  double unitScale = 138.935485;
  double direrror, recerror, tmperror;
  ele.reinit (beta, K, rcut, *box, method, order);
  ele.errorEstimatePME (c2, N, direrror, tmperror);
  ele.errorEstimateInterpolRec (c2, N, recerror);
  double dirtime, looptime, convtime;
  tc.cal_time (rcut, K, order, dirtime, looptime, convtime);
  std::cout << "**************************************************\n"
	    << "** beta " << beta
	    << "\n** rcut " << rcut
	    << "\n** K " << K[0]  << " " << K[1] << " " << K[2]
	    << "\n** order " << order
	    << "\n** required error is " << higherRequire * 2 * unitScale
	    << "\tdirerror " << direrror * unitScale
	    << "\trecerror " << recerror * unitScale
	    << "\tsum is " << (direrror + recerror) * unitScale
	    << "\n** time dir " << dirtime 
	    << "\tloop " << looptime 
	    << "\tconv " << convtime 
	    << "\trec " << looptime + convtime 
	    << "\ttotal " << dirtime + looptime + convtime 
	    << "\n**************************************************"
	    << std::endl;
  return true;
}


// template <typename Judge >
// double NewtonMin (const double & start,
// 		  Judge & f, 
// 		  const double & precision)
// {
//   double next = start;
//   double now = start - 2 * precision;
//   double h = 0.001;
//   double valueNow0, valueNow1, valueNow2;
//   int count  = 0;
  
//   while (fabs(now - next) > precision && count ++ < 10){
    
//     std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
//     std::cout << "%%% step " << count 
// 	      << "\tnow at " << now  << "\tnext is " << next << std::endl;
//     std::cout << "%%% newton Min, the diff step is " << fabs(now - next) << std::endl;
//     now = next;
//     valueNow0 = f(now - h);
//     std::cout << "rcut " << f.rcut
// 	      << "\torder " << f.order
// 	      << "\tK " << f.K[0]  << " " << f.K[1] << " " << f.K[2]
// 	      << "\tdir " << f.dirtime 
// 	      << "\tloop " << f.looptime 
// 	      << "\tconv " << f.convtime 
// 	      << "\trec " << f.looptime + f.convtime 
// 	      << "\ttotal " << f.dirtime + f.looptime + f.convtime 
// 	      << std::endl;
//     valueNow1 = f(now);
//     std::cout << "rcut " << f.rcut
// 	      << "\torder " << f.order
// 	      << "\tK " << f.K[0]  << " " << f.K[1] << " " << f.K[2]
// 	      << "\tdir " << f.dirtime 
// 	      << "\tloop " << f.looptime 
// 	      << "\tconv " << f.convtime 
// 	      << "\trec " << f.looptime + f.convtime 
// 	      << "\ttotal " << f.dirtime + f.looptime + f.convtime 
// 	      << std::endl;
//     valueNow2 = f(now + h);
//     std::cout << "rcut " << f.rcut
// 	      << "\torder " << f.order
// 	      << "\tK " << f.K[0]  << " " << f.K[1] << " " << f.K[2]
// 	      << "\tdir " << f.dirtime 
// 	      << "\tloop " << f.looptime 
// 	      << "\tconv " << f.convtime 
// 	      << "\trec " << f.looptime + f.convtime 
// 	      << "\ttotal " << f.dirtime + f.looptime + f.convtime 
// 	      << std::endl;
//     double grad = (valueNow2 - valueNow0) / (2*h);
//     double hess = (valueNow2 - 2 * valueNow1 + valueNow0) / h / h;
//     next = now - grad / hess;
//     if (now < 0){
// //       now = 0.5;
//     }
//   }
//   return f(next);
// }
    
