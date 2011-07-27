#ifndef __AssignForceCorr_h_wanghan__
#define __AssignForceCorr_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDSystem_interface.h"
#include "ForceCorr.h"

class AssignForceCorr 
{
  RectangularBox box;
  int  nx, ny, nz, nele;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  ScalorType * hfcx, * hfcy, *hfcz;
  ScalorType * dfcx, * dfcy, *dfcz;
  bool malloced;
  dim3 myBlockDim;
  dim3 atomGridDim;
  mutable FILE * fp_write;
private:
  void freeAll ();
public:
  AssignForceCorr ();
  ~AssignForceCorr ();
  void reinit (const MDSystem & sys,
	       const ForceCorr & arc,
	       const IndexType & NThread);
  void getForceCorr (const ForceCorr & arc);
  void zero ();
  void assign (MDSystem & sys);
public:
  void print_x  (const char * file) const;
  void init_write (const char * file) const;
  void write (const ScalorType & time) const;
  void end_write () const;
}
    ;

unsigned AssignForceCorr::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void AssignForceCorr::
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
