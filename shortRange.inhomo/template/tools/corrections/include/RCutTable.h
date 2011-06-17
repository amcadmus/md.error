#ifndef __RCutTable_h_wanghan__
#define __RCutTable_h_wanghan__

#include "DensityProfile.h"

class RCutTable
{
  std::vector<double > boxsize;
  int  nx, ny, nz, nele;
  double hx, hy, hz;
  int nrc;
  std::vector<double > rcList;
  std::vector<int    > rcutIndex;
public:
  RCutTable ();
  ~RCutTable ();
  void load_rc (const char * file);
  void save_rc (const char * file) const;
  const std::vector<int    > & getProfileIndex () const {return rcutIndex;}
  const std::vector<double > & getRcList       () const {return rcList;}
  std::vector<int > & getProfileIndex () {return rcutIndex;}
};


#endif
