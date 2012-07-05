#include "Stopwatch.h"
#include <iostream>

double Stopwatch::user () const
{
  return (double)(toc.tms_utime - tic.tms_utime) * HZi;
}

double Stopwatch::system () const
{
  return (double)(toc.tms_stime - tic.tms_stime) * HZi;
}

double Stopwatch::real () const
{
  return (double)(r1 - r0) * HZi;
}

void Stopwatch::stop ()
{
  r1 = times (&toc);
}

void Stopwatch::start() 
{
  r0 = times (&tic);
}

Stopwatch::Stopwatch ()
    : HZi (1./HZ)
{
}
