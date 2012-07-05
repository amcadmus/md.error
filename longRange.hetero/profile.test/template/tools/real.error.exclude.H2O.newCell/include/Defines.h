#ifndef __Defines_h_wanghan__
#define __Defines_h_wanghan__

typedef double ValueType ;

struct IntVectorType
{
  int x, y, z;
  IntVectorType (const int xx=0, const int yy=0, const int zz=0);
};

struct VectorType
{
  double x, y, z;
  VectorType (const double xx=0, const double yy=0, const double zz=0);
};

inline IntVectorType::
IntVectorType (const int xx, const int yy, const int zz)
    : x(xx), y(yy), z(zz)
{
}

inline VectorType::
VectorType (const double xx, const double yy, const double zz)
    : x(xx), y(yy), z(zz)
{
}

#endif
