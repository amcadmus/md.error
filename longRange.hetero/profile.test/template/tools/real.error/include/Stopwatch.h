#ifndef __Stopwatch_h_wanghan__
#define __Stopwatch_h_wanghan__

#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>
#ifndef HZ
#define HZ 100.0
#endif

class Stopwatch 
{ 
public: 
  Stopwatch(); 
 
  void              start(); 
  void              stop(); 
  double            system() const; 
  double            user() const; 
  double            real() const; 
 
  static double     resolution() {return 1./HZ;}; 
private:
  struct tms tic, toc;
  long r1, r0;
  double HZi;
};               


#endif
// end of file 

