#ifndef __AssignRCut_h_wanghan__
#define __AssignRCut_h_wanghan__

#define DEVICE_CODE

#include "common.h"
#include "MDSystem_interface.h"

class AssignRCut 
{
  RectangularBox box;
  int  nx, ny, nz, nele;
  // double hx, hy, hz;
  ScalorType * hrcut;
  ScalorType * drcut;
  bool malloced;
  dim3 myBlockDim;
  dim3 atomGridDim;
private:
  void freeAll ();
public:
  AssignRCut ();
  ~AssignRCut ();
  void reinit (const char * filename,
	       const MDSystem & sys,
	       const IndexType & NThread);
  void assign (MDSystem & sys);
}
    ;

#endif
