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

#include "GroFileManager.h"

using namespace std;

int main(int argc, char * argv[])
{
  string ifile ("conf.gro");
  string ofile ("out.gro");
  vector<int > resdindex, atomindex;
  vector<string > resdname, atomname;
  vector<vector<double > > posi, velo;
  vector<double > boxsize;


  GroFileManager::read (ifile, resdindex, resdname, atomname, atomindex,
			posi, velo, boxsize);
  for (unsigned i = 0; i < posi.size(); ++i){
    // cout << i << endl;
    posi[i][0] += .5 * boxsize[0];
  }
  boxsize[0] *= 2.;
  GroFileManager::write (ofile, resdindex, resdname, atomname, atomindex,
			 posi, velo, boxsize);
  
  return 0;
}
