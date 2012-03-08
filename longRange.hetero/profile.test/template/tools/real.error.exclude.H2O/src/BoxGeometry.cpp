#include "BoxGeometry.h"

TimeConstRectangularBoxGeometry::TimeConstRectangularBoxGeometry (const double & l0, 
								  const double & l1,
								  const double & l2)
{
  boxl.clear();
  li.clear();
  boxl.resize (3, 0);
  li.resize (3, 0);
  boxl[0] = l0;
  boxl[1] = l1;
  boxl[2] = l2;
  li[0] = 1./l0;
  li[1] = 1./l1;
  li[2] = 1./l2;
}

void TimeConstRectangularBoxGeometry::getBaseVector (std::vector<std::vector<value_type > > & baseVector)
    const
{
  baseVector.clear();
  std::vector<double > tmp(3, 0);
  tmp[0] = boxl[0];
  baseVector.push_back (tmp);
  tmp[0] = 0;
  tmp[1] = boxl[1];
  baseVector.push_back (tmp);
  tmp[1] = 0;
  tmp[2] = boxl[2];
  baseVector.push_back (tmp);
}

void TimeConstRectangularBoxGeometry::getBoxSize (std::vector<value_type > & boxl_) 
    const
{
  boxl_ = boxl;
}


inline//  void TimeConstRectangularBoxGeometry::moveParticleToBox (const double & time,
// 						       StandardParticle & p) const
// {
//   if (p.r()[0] >= 0){
//     p.noi()[0] -= int(p.r()[0] * li[0] + .5);
//   }
//   else{
//     p.noi()[0] -= int(p.r()[0] * li[0] - .5);
//   }
//   p.r()[0] += p.noi()[0] * boxl[0];

//   if (p.r()[1] >= 0){
//     p.noi()[1] -= int(p.r()[1] * li[1] + .5);
//   }
//   else{
//     p.noi()[1] -= int(p.r()[1] * li[1] - .5);
//   }
//   p.r()[1] += p.noi()[1] * boxl[1];

//   if (p.r()[2] >= 0){
//     p.noi()[2] -= int(p.r()[2] * li[2] + .5);
//   }
//   else{
//     p.noi()[2] -= int(p.r()[2] * li[2] - .5);
//   }
//   p.r()[2] += p.noi()[2] * boxl[2];
// }

void TimeConstRectangularBoxGeometry::moveParticleToBox (std::vector<StandardParticle > & p) const
{
  for (std::vector<StandardParticle >::iterator q = p.begin(); q != p.end(); q ++){
    int tmp;
    if ((*q).r()[0] >= 0){
      tmp = int((*q).r()[0] * li[0] + .5);
      (*q).noi()[0] -= tmp;
    }
    else{
      tmp = int((*q).r()[0] * li[0] - .5);
      (*q).noi()[0] -= tmp;
    }
    (*q).r()[0] -= tmp * boxl[0];

    if ((*q).r()[1] >= 0){
      tmp = int((*q).r()[1] * li[1] + .5);
      (*q).noi()[1] -= tmp;
    }
    else{
      tmp = int((*q).r()[1] * li[1] - .5);
      (*q).noi()[1] -= tmp;
    }
    (*q).r()[1] -= tmp * boxl[1];

    if ((*q).r()[2] >= 0){
      tmp = int((*q).r()[2] * li[2] + .5);
      (*q).noi()[2] -= tmp;
    }
    else{
      tmp = int((*q).r()[2] * li[2] - .5);
      (*q).noi()[2] -= tmp;
    }
    (*q).r()[2] -= tmp * boxl[2];
  }
}


// inline void TimeConstRectangularBoxGeometry::shortestImage (const double & time, 
// 						   std::vector<value_type > & diff) const 
// {
//   int i = 0;
//   if (diff[i] > 0.5 * boxl[i]){
//     diff[i] -= boxl[i];
//   }
//   else if (diff[i] < - 0.5 * boxl[i]){
//     diff[i] += boxl[i];
//   }
//   i = 1;
//   if (diff[i] > 0.5 * boxl[i]){
//     diff[i] -= boxl[i];
//   }
//   else if (diff[i] < - 0.5 * boxl[i]){
//     diff[i] += boxl[i];
//   }
//   i = 2;
//   if (diff[i] > 0.5 * boxl[i]){
//     diff[i] -= boxl[i];
//   }
//   else if (diff[i] < - 0.5 * boxl[i]){
//     diff[i] += boxl[i];
//   }
// }

    
void RectangularBoxGeometry::setBoxSize (const std::vector<value_type > & boxl_)
{
  boxl = boxl_;
}

void RectangularBoxGeometry::setBoxSize (const double & a, 
					 const double & b, 
					 const double & c) 
{
  boxl.clear();
  boxl.push_back (a);
  boxl.push_back (b);
  boxl.push_back (c);
}


void RectangularBoxGeometry::getBaseVector (std::vector<std::vector<value_type > > & baseVector)
    const
{
  baseVector.clear();
  std::vector<double > tmp(0, 3);
  tmp[0] = boxl[0];
  baseVector.push_back (tmp);
  tmp[0] = 0;
  tmp[1] = boxl[1];
  baseVector.push_back (tmp);
  tmp[1] = 0;
  tmp[2] = boxl[2];
  baseVector.push_back (tmp);
}

void RectangularBoxGeometry::getBoxSize (std::vector<value_type > & boxl_) 
    const
{
  boxl_ = boxl;
}
