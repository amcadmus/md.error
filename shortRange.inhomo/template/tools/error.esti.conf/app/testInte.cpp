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

#include "Integral3D.h"
#include "Integral1D.h"
#include "Test.h"

int main(int argc, char * argv[])
{
  calculateGlobal ();
  double rc = 5.;
  double prec = 1e-11;  

  InteNaive00 fi00;
  Integral3D<InteNaive00, double > i00;
  double value_i00 = i00.cal_int(Integral3DInfo::Gauss27,
				 fi00,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive11 fi11;
  Integral3D<InteNaive11, double > i11;
  double value_i11 = i11.cal_int(Integral3DInfo::Gauss27,
				 fi11,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive22 fi22;
  Integral3D<InteNaive22, double > i22;
  double value_i22 = i22.cal_int(Integral3DInfo::Gauss27,
				 fi22,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive01 fi01;
  Integral3D<InteNaive01, double > i01;
  double value_i01 = i01.cal_int(Integral3DInfo::Gauss27,
				 fi01,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive02 fi02;
  Integral3D<InteNaive02, double > i02;
  double value_i02 = i02.cal_int(Integral3DInfo::Gauss27,
				 fi02,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive12 fi12;
  Integral3D<InteNaive12, double > i12;
  double value_i12 = i12.cal_int(Integral3DInfo::Gauss27,
				 fi12,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);

  InteNaive0 fi0;
  Integral3D<InteNaive0, double > i0;
  double value_i0 = i0.cal_int(Integral3DInfo::Gauss27,
				 fi0,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive1 fi1;
  Integral3D<InteNaive1, double > i1;
  double value_i1 = i1.cal_int(Integral3DInfo::Gauss27,
				 fi1,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);
  InteNaive2 fi2;
  Integral3D<InteNaive2, double > i2;
  double value_i2 = i2.cal_int(Integral3DInfo::Gauss27,
				 fi2,
				 rc, 50,
				 0, 2*M_PI,
				 0, M_PI,
				 prec);

  double b0k, b1k, b2k;
  double a00k, a01k, a02k, a11k, a12k, a22k;
  inteClever (1, 1, rc,
	      b0k, b1k, b2k,
	      a00k, a01k, a02k, a11k, a12k, a22k);
  
  printf ("i00_3d: %e, i00_1d: %e, diff: %e\n",
	  value_i00, a00k, fabs(a00k - value_i00));
  printf ("i11_3d: %e, i11_1d: %e, diff: %e\n",
	  value_i11, a11k, fabs(a11k - value_i11));
  printf ("i22_3d: %e, i22_1d: %e, diff: %e\n",
	  value_i22, a22k, fabs(a22k - value_i22));
  printf ("i01_3d: %e, i01_1d: %e, diff: %e\n",
	  value_i01, a01k, fabs(a01k - value_i01));
  printf ("i02_3d: %e, i02_1d: %e, diff: %e\n",
	  value_i02, a02k, fabs(a02k - value_i02));
  printf ("i12_3d: %e, i12_1d: %e, diff: %e\n",
	  value_i12, a12k, fabs(a12k - value_i12));
  printf ("i0_3d: %e, i0_1d: %e, diff: %e\n",
	  value_i0, b0k, fabs(b0k - value_i0));
  printf ("i1_3d: %e, i1_1d: %e, diff: %e\n",
	  value_i1, b1k, fabs(b1k - value_i1));
  printf ("i2_3d: %e, i2_1d: %e, diff: %e\n",
	  value_i2, b2k, fabs(b2k - value_i2));
			       
}

