#ifndef __wanghan_Particle_h__
#define __wanghan_Particle_h__

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>

typedef double value_type;
typedef long int id_type;
typedef long int type_type;
#define NDIM 3

class PositionParticle 
{
public:
  virtual ~PositionParticle () {};
  virtual std::vector<value_type > & r (void) = 0;
  virtual const std::vector<value_type > & r (void) const = 0;
  virtual id_type & id (void) = 0;
  virtual const id_type & id (void) const = 0;
};

class StandardParticle : public PositionParticle {
  std::vector<value_type > rBuff;	// position
  std::vector<value_type > vBuff;	// velocity
  std::vector<value_type > fBuff;	// force
  std::vector<value_type > noiBuff;	// number of images
  value_type massBuff;			// mass
  value_type chargeBuff;		// charge
  id_type idBuff;			// id
  type_type typeBuff;			// type
public:
  StandardParticle ();
  StandardParticle (const StandardParticle & part);
  ~StandardParticle () {};
public:
  void clear ();
  std::vector<value_type > & r (void) {return rBuff;}
  const std::vector<value_type > & r (void) const {return rBuff;}
  std::vector<value_type > & v (void) {return vBuff;}
  const std::vector<value_type > & v (void) const {return vBuff;}
  std::vector<value_type > & f (void) {return fBuff;}
  const std::vector<value_type > & f (void) const {return fBuff;}
  std::vector<value_type > & noi (void) {return noiBuff;}
  const std::vector<value_type > & noi (void) const {return noiBuff;}
  value_type & mass (void) {return massBuff;}
  const value_type & mass (void) const {return massBuff;}
  value_type & charge (void) {return chargeBuff;}
  const value_type & charge (void) const {return chargeBuff;}
  id_type & id (void) {return idBuff;}
  const id_type & id (void) const {return idBuff;}
  type_type & type (void) {return typeBuff;}
  const type_type & type (void) const {return typeBuff;}
};


#endif
