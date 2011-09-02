#ifndef __AssignFcorr_h_wanghan__
#define __AssignFcorr_h_wanghan__

#include <stdlib.h>
#include <iostream>


class AssignFcorr 
{
  double box[3];
  int  nx, ny, nz, nele;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  float *hfx, *hfy, *hfz;
  bool malloced;
  mutable FILE * fp_write;
  mutable FILE * fp_read;
private:
  void freeAll ();
public:
  AssignFcorr ();
  ~AssignFcorr ();
  void uniform (const double & rc);
public:
  void init_read	 (const char * file);
  bool read		 (float & time);
  void end_read		 () const;
public:
  void print_xy		 (const char * file) const;
}
    ;

unsigned AssignFcorr::
index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const
{
  return iz + nz * (iy + ny * ix);
}

void AssignFcorr::
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
