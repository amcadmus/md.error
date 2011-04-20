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

#include "ErrorProfile.h"

int main(int argc, char * argv[])
{
  if (argc != 3){
    std::cerr << "usage:\n"
	     << argv[0] << " forceFile refh" << std::endl;
    return 1;
  }

  ErrorProfile_PiecewiseConst dp (std::string(argv[1]), atof(argv[2]));

  dp.print_x (std::string("real.x.out"));
  dp.print_xy(std::string("real.xy.out"));
  
  return 0;
}


