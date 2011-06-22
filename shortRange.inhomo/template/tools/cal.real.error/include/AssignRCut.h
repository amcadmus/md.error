#ifndef __AssignRCut_h_wanghan__
#define __AssignRCut_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDSystem_interface.h"
#include "AdaptRCut.h"

class AssignRCut 
{
  RectangularBox box;
  int  nx, ny, nz, nele;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  ScalorType * hrcut;
  ScalorType * drcut;
  bool malloced;
  dim3 myBlockDim;
  dim3 atomGridDim;
  mutable FILE * fp_write;
private:
  void freeAll ();
public:
  AssignRCut ();
  ~AssignRCut ();
  void reinit (const MDSystem & sys,
	       const AdaptRCut & arc,
	       const IndexType & NThread);
  void getRCut (const AdaptRCut & arc);
  void uniform (const double & rc);
  void assign (MDSystem & sys);
public:
  void print_x  (const char * file) const;
  void init_write (const char * file) const;
  void write (const ScalorType & time) const;
  void end_write () const;
}
    ;

unsigned AssignRCut::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void AssignRCut::
index1to3 (unsigned& input,
	   unsigned& ix, unsigned& iy, unsigned& iz) const
{
  unsigned tmp = input;
  iz = tmp % (nz);
  tmp = (tmp - iz) / nz;
  iy = tmp % (ny);
  ix =  (tmp - iy) / ny;
}

#endif
