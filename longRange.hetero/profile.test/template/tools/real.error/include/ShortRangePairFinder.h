#ifndef __wanghan_ShortRangePairFinder_h__
#define __wanghan_ShortRangePairFinder_h__

#include "Particle.h"
#include "BoxGeometry.h"
#include <vector>
#include <algorithm>
#include <list>

typedef PositionParticle * pPositionParticle;
typedef double value_type;

class ShortRangePairFinder 
{
public:
  typedef std::vector<PositionParticle * >::iterator		PairListIterator;
  typedef std::vector<PositionParticle * >::const_iterator	const_PairListIterator;
public:
  virtual ~ShortRangePairFinder () {};
  virtual PairListIterator 		begin() = 0;
  virtual const_PairListIterator 	begin() const = 0;
  virtual PairListIterator		end() = 0;
  virtual const_PairListIterator 	end() const = 0;
}
    ;


class NeighborListPairFinder : public ShortRangePairFinder
{
private:
  value_type rcut;
  value_type rshell;
  value_type rtotal;
  value_type rrtotal;
  // these are setted during construction
  const TimeConstRectangularBoxGeometry * pbox;
  std::vector<value_type > boxl;
  std::vector<unsigned > nCell;
  std::vector<value_type > celll;
  std::vector<value_type > cellli;
  std::vector<std::vector<int > > switchVecs;
  // these will vary during simulation
  std::vector<std::vector<std::vector<value_type > > > oldPositions;
  std::vector<std::list<PositionParticle * > > lists;
  std::vector<PositionParticle * > neighborList;
private:
  bool checkMove ();
public:
  NeighborListPairFinder (const double rcutout, const double rshellout,
			  const TimeConstRectangularBoxGeometry & box);
  // usr should garentee that the particle is in the box before
  // registration
  void registerParticle (PositionParticle & p);
  void registerParticle (std::vector<PositionParticle * > & p);
  void unregisterParticle (PositionParticle & p);
  void unregisterParticle (std::vector<PositionParticle * > & p);
  // if all registered particles are in the right cell, then use
  // build. generally speaking, use build() at first time when
  // particles are registered, use rebuild() if some particles are
  // unregistered and use checkAndRebuild() during simulation.
  void build ();
  // if some registered particles move out of their cell, then use
  // rebuild (we rearrange out-of-cell particles first)
  void rebuild ();
  // check whether the neighbor list is valid and if not rebuild the
  // list.
  void checkList ();
public:
  PairListIterator begin () {return neighborList.begin();}
  const_PairListIterator begin () const {return neighborList.begin();}
  PairListIterator end () {return neighborList.end();}
  const_PairListIterator end () const {return neighborList.end();}
}
    ;

#endif
 
//       {	
// 	if (iter1 != list.end()){
// 	  p0 = (*iter0);
// 	  p1 = (*iter1);
// 	  iter1 ++;
// 	  return true;
// 	}
// 	else {
// 	  iter0 ++;
// 	  preiter0 ++;
// 	  if (preiter0 == list.end()){
// 	    return false;
// 	  }
// 	  else {
// 	    iter1 = preiter0;
// 	    p0 = (*iter0);
// 	    p1 = (*iter1);
// 	    iter1 ++;
// 	    return true;
// 	  }
// 	}
//       }
