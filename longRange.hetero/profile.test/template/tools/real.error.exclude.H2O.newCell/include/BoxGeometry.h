#ifndef __wanghan_BoxGeometry_h__
#define __wanghan_BoxGeometry_h__

#include <vector>
#include "Particle.h"

typedef double value_type;

class BoxGeometry
{
public:
  virtual ~BoxGeometry () {};
  virtual void moveParticleToBox (StandardParticle & p) const=0;
  virtual void moveParticleToBox (std::vector<StandardParticle > & p) const=0;
  virtual void shortestImage (std::vector<value_type> & diff) const=0;
}
    ;

class TriclinicBoxGeometry : public BoxGeometry
{
public:
  virtual ~TriclinicBoxGeometry () {};
  virtual void getBaseVector (std::vector<std::vector<value_type > > & baseVector)
      const = 0;
};

class RectangularBoxGeometry : public TriclinicBoxGeometry
{
protected:
  std::vector<value_type > boxl;
  std::vector<value_type > li;
public:
  RectangularBoxGeometry () {};
  ~RectangularBoxGeometry () {};
public:
  void setBoxSize (const std::vector<value_type > & boxl);
  void setBoxSize (const double & a, 
		   const double & b, 
		   const double & c);
public:
  void getBaseVector(std::vector<std::vector<value_type > > & baseVector)
      const;
  void getBoxSize (std::vector<value_type > & boxl) 
      const;
};


class TimeConstRectangularBoxGeometry : public TriclinicBoxGeometry
{
private:
  std::vector<value_type > boxl;
  std::vector<value_type > li;
public:
  TimeConstRectangularBoxGeometry (const double & l0, 
				   const double & l1,
				   const double & l2);
  ~TimeConstRectangularBoxGeometry () {};
public:
  void getBaseVector(std::vector<std::vector<value_type > > & baseVector)
      const;
  void getBoxSize (std::vector<value_type > & boxl) 
      const;
public:
  inline void moveParticleToBox (StandardParticle & p) const{
    int tmp;
    if (p.r()[0] >= 0){
      tmp = int(p.r()[0] * li[0] + .5);
      p.noi()[0] -= tmp;
    }
    else{
      tmp = int(p.r()[0] * li[0] - .5);
      p.noi()[0] -= tmp;
    }
    p.r()[0] -= tmp * boxl[0];

    if (p.r()[1] >= 0){
      tmp = int(p.r()[1] * li[1] + .5);
      p.noi()[1] -= tmp;
    }
    else{
      tmp = int(p.r()[1] * li[1] - .5);
      p.noi()[1] -= tmp;
    }
    p.r()[1] -= tmp * boxl[1];

    if (p.r()[2] >= 0){
      tmp = int(p.r()[2] * li[2] + .5);
      p.noi()[2] -= tmp;
    }
    else{
      tmp = int(p.r()[2] * li[2] - .5);
      p.noi()[2] -= tmp;
    }
    p.r()[2] -= tmp * boxl[2];
  }
  void moveParticleToBox (std::vector<StandardParticle > & p) const;
  inline void shortestImage (std::vector<value_type> & diff) const {
    int i = 0;
    if (diff[i] >= 0.5 * boxl[i]){
      diff[i] -= boxl[i];
    }
    else if (diff[i] < - 0.5 * boxl[i]){
      diff[i] += boxl[i];
    }
    i = 1;
    if (diff[i] >= 0.5 * boxl[i]){
      diff[i] -= boxl[i];
    }
    else if (diff[i] < - 0.5 * boxl[i]){
      diff[i] += boxl[i];
    }
    i = 2;
    if (diff[i] >= 0.5 * boxl[i]){
      diff[i] -= boxl[i];
    }
    else if (diff[i] < - 0.5 * boxl[i]){
      diff[i] += boxl[i];
    }
  }
}
    ;


#endif
