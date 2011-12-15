#ifndef __AssignRCut_h_wanghan__
#define __AssignRCut_h_wanghan__

#include <stdlib.h>
#include <iostream>


class AssignRCut 
{
  double box[3];
  int  nx, ny, nz, nele;
  inline unsigned index3to1 (unsigned  ix, unsigned  iy, unsigned  iz) const;
  inline void     index1to3 (unsigned& input,
			     unsigned& ix, unsigned& iy, unsigned& iz) const;
  float * hrcut;
  bool malloced;
  mutable FILE * fp_write;
  mutable FILE * fp_read;
private:
  void freeAll ();
public:
  AssignRCut ();
  ~AssignRCut ();
  void uniform (const double & rc);
public:
  void init_write	 (const char * file) const;
  void write		 (const float & time) const;
  void end_write	 () const;
public:
  void init_read	 (const char * file);
  bool read		 (float & time);
  void end_read		 () const;
public:
  void print_x		 (const char * file) const;
  void print_xy		 (const char * file) const;
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
