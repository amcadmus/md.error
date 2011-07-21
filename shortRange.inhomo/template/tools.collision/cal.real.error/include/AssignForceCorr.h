#ifndef __AssignForceCorr_h_wanghan__
#define __AssignForceCorr_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDSystem_interface.h"

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
  mutable FILE * fp_read;
private:
  void freeAll ();
public:
  AssignForceCorr ();
  ~AssignForceCorr ();
  void reinit (const MDSystem & sys,
	       const IndexType & NThread);
  void zero ();
  inline void getForceCorr (const ScalorType & x,
			    const ScalorType & y,
			    const ScalorType & z,
			    ScalorType & fcx,
			    ScalorType & fcy,
			    ScalorType & fcz) const;
public:
  void print_x  (const char * file) const;
  void init_write (const char * file) const;
  void write (const ScalorType & time) const;
  void end_write () const;
  void init_read (const char * file);
  void read (ScalorType & time);
  void end_read () const;
}
    ;

void AssignForceCorr::
getForceCorr (const ScalorType & x,
	      const ScalorType & y,
	      const ScalorType & z,
	      ScalorType & fcx,
	      ScalorType & fcy,
	      ScalorType & fcz) const
{
  ScalorType tmpx (x), tmpy(y), tmpz(z);
  if      (tmpx <  0         ) tmpx += box.size.x;
  else if (tmpx >= box.size.x) tmpx -= box.size.x;
  if      (tmpy <  0         ) tmpy += box.size.y;
  else if (tmpy >= box.size.y) tmpy -= box.size.y;
  if      (tmpz <  0         ) tmpz += box.size.z;
  else if (tmpz >= box.size.z) tmpz -= box.size.z;
  IndexType idx = index3to1(tmpx / box.size.x * nx,
			    tmpy / box.size.y * ny ,
			    tmpz / box.size.z * nz);
  fcx = hfcx[idx];
  fcy = hfcy[idx];
  fcz = hfcz[idx];
}


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
