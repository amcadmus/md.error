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
  Integral3D<tmpf3, double > i3;
  tmpf3 f3;
  double prec = 1e-6;
  
  double value_i3 = i3.cal_int(Integral3DInfo::Gauss27,
			       f3,
			       5, 20,
			       0, 2*M_PI,
			       0, M_PI,
			       prec);

  Integral1D<tmpf, double > i1;
  tmpf f1;
  double value_i1 = i1.cal_int(Integral1DInfo::Gauss4,
			       f1,
			       5, 20,
			       prec);
  
  printf ("i3: %e, i1: %e, diff: %e\n",
	  value_i3, value_i1, fabs(value_i1 - value_i3));
			       
}

