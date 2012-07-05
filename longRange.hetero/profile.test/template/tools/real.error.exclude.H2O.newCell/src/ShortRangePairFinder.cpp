#include "ShortRangePairFinder.h"
#include "VectorOperation.h"

void makeSwitchVector (std::vector<std::vector<int > > & vecs);

NeighborListPairFinder::NeighborListPairFinder (const double rcutout, 
						const double rshellout,
						const TimeConstRectangularBoxGeometry & box)
    : rcut (rcutout),
      rshell (rshellout),
      rtotal (rcut+rshell),
      rrtotal (rtotal*rtotal),
      pbox (&box),
      boxl (3, 0),
      nCell (3, 0),
      celll (3, 0),
      cellli (3, 0)
{
  pbox->getBoxSize(boxl);
  nCell[0] = unsigned (boxl[0] / (rtotal));
  nCell[1] = unsigned (boxl[1] / (rtotal));
  nCell[2] = unsigned (boxl[2] / (rtotal));
  celll[0] = boxl[0] / nCell[0];
  celll[1] = boxl[1] / nCell[1];
  celll[2] = boxl[2] / nCell[2];
  cellli[0] = 1./celll[0];
  cellli[1] = 1./celll[1];
  cellli[2] = 1./celll[2];
  lists.resize ((nCell[0]) * (nCell[1]) * (nCell[2]));
  makeSwitchVector (switchVecs);
}

// static int count = 0;
void NeighborListPairFinder::registerParticle (PositionParticle & p)
{
  std::vector<unsigned > index(3);
  if (p.r()[0] + 0.5 * boxl[0] == boxl[0]) {p.r()[0] = - 0.5 * boxl[0]; /*std::cout << "" << std::endl;*/}
  if (p.r()[1] + 0.5 * boxl[1] == boxl[1]) {p.r()[1] = - 0.5 * boxl[1]; /*std::cout << "" << std::endl;*/}
  if (p.r()[2] + 0.5 * boxl[2] == boxl[2]) {p.r()[2] = - 0.5 * boxl[2]; /*std::cout << "" << std::endl;*/}
  index[0] = unsigned((p.r()[0] + 0.5 * boxl[0]) * cellli[0]);
  index[1] = unsigned((p.r()[1] + 0.5 * boxl[1]) * cellli[1]);
  index[2] = unsigned((p.r()[2] + 0.5 * boxl[2]) * cellli[2]);
//  std::cout << count ++ << "\t" << lists[index[0] + index[1]*nCell[0] + index[2]*nCell[0]*nCell[1]].size() << std::endl;
  lists[index[0] + index[1]*nCell[0] + index[2]*nCell[0]*nCell[1]].push_back(&p);
}


void NeighborListPairFinder::registerParticle (std::vector<PositionParticle * > & p)
{
  std::vector<unsigned > index(3);
  for (std::vector<PositionParticle * >::iterator pp = p.begin(); 
       pp != p.end(); pp ++){    
    index[0] = unsigned(((*pp)->r()[0] + 0.5 * boxl[0]) * cellli[0]);
    index[1] = unsigned(((*pp)->r()[1] + 0.5 * boxl[1]) * cellli[1]);
    index[2] = unsigned(((*pp)->r()[2] + 0.5 * boxl[2]) * cellli[2]);
    lists[index[0] + index[1]*nCell[0] + index[2]*nCell[0]*nCell[1]].push_back(&(*(*pp)));
  }
}

void NeighborListPairFinder::rebuild (void)
{  
  // first rearrange the particles
  for (unsigned k = 0; k < nCell[2]; k ++){
    for (unsigned j = 0; j < nCell[1]; j ++){
      for (unsigned i = 0; i < nCell[0]; i ++){
	unsigned thisIndex = i + j * nCell[0] + k * nCell[0] * nCell[1];
	// iter over this list
	std::vector<unsigned > partCood (3);
	for (std::list<PositionParticle * >::iterator it = lists[thisIndex].begin();
	     it != lists[thisIndex].end(); ){
	  partCood[0] = unsigned(((*it)->r()[0] + 0.5 * boxl[0]) * cellli[0]);
	  partCood[1] = unsigned(((*it)->r()[1] + 0.5 * boxl[1]) * cellli[1]);
	  partCood[2] = unsigned(((*it)->r()[2] + 0.5 * boxl[2]) * cellli[2]);
	  unsigned partIndex = partCood[0] + partCood[1]*nCell[0] + partCood[2]*nCell[0]*nCell[1];
	  if (partIndex != thisIndex){
	    std::list<PositionParticle * >::iterator del = it;
	    ++ it;
	    lists[partIndex].push_back (*del);
	    lists[thisIndex].erase (del);
	  }
	  else{
	    ++ it;
	  }
	}
      }
    }
  }
  // then build as initially
  build();
}

	    

void NeighborListPairFinder::build (void)
{
  // clear neighborList
  neighborList.clear();
  // remember old positions
  oldPositions.clear();
  oldPositions.resize((nCell[0]) * (nCell[1]) * (nCell[2]));
  std::vector<std::list<PositionParticle * > >::iterator plist = lists.begin();
  std::vector<std::vector<std::vector<value_type > > >::iterator ppos = oldPositions.begin();
//   std::cout << lists.size() << std::endl;
//   std::cout << oldPositions.size()  << std::endl;
//   int count = 0;
  for (; plist != lists.end(); plist++, ppos ++){
    for (std::list<PositionParticle * >::iterator ppart = (*plist).begin();
	 ppart != (*plist).end(); ppart ++){
      (*ppos).push_back ((*ppart)->r());
    }
//     std::cout << count ++ << '\t' << ppos->size() << std::endl;
  }
//   std::cout << std::endl;

  for (unsigned i = 0; i < nCell[0]; i ++){
    for (unsigned j = 0; j < nCell[1]; j ++){
      for (unsigned k = 0; k < nCell[2]; k ++){
	for (std::vector<std::vector<int > >::iterator sw = switchVecs.begin();
	     sw != switchVecs.end(); sw ++){
	  unsigned thisCell = i + j * nCell[0] + k * nCell[0] * nCell[1];
	  int tmp;
	  unsigned nextCell = 0;
	  std::vector<value_type > nextCellShift (3);
	  tmp = i + (*sw)[0];
	  if (tmp ==-1){
	    tmp = nCell[0]-1;
	    nextCellShift[0] = -boxl[0];
	  }
	  else if (tmp == int(nCell[0])){
	    tmp = 0;
	    nextCellShift[0] = boxl[0];
	  }
	  nextCell += tmp;
	  tmp = j + (*sw)[1];
	  if (tmp == -1){
	    tmp = nCell[1]-1;
	    nextCellShift[1] = -boxl[1];
	  }
	  else if (tmp == int(nCell[1])){
	    tmp = 0;
	    nextCellShift[1] = boxl[1];
	  }
	  nextCell += tmp * nCell[0];
	  tmp = k + (*sw)[2];
	  if (tmp == -1){
	    tmp = nCell[2]-1;
	    nextCellShift[2] = -boxl[2];
	  }
	  else if (tmp == int(nCell[2])){
	    tmp = 0;
	    nextCellShift[2] = boxl[2];
	  }
	  nextCell += tmp * nCell[0] * nCell[1];
	  
	  
	  for (std::list<PositionParticle * >::iterator p0 = lists[thisCell].begin();
	       p0 != lists[thisCell].end(); p0 ++){
	    for (std::list<PositionParticle *>::iterator p1 = lists[nextCell].begin();
		 p1 != lists[nextCell].end(); p1 ++){
	      if (thisCell != nextCell || (*p0)->id() < (*p1)->id()){
		std::vector<double > diff ((*p0)->r());
		VectorOperation::add (diff, -1, (*p1)->r(), -1, nextCellShift);
		if (VectorOperation::dot (diff, diff) < rrtotal){
		  neighborList.push_back(*p0);
		  neighborList.push_back(*p1);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

bool NeighborListPairFinder::checkMove ()
{
  value_type value = 0.25 * rshell * rshell;
  std::vector<std::list<PositionParticle * > >::iterator itlists = lists.begin();
  std::vector<std::vector<std::vector<value_type > > >::iterator itposlists = oldPositions.begin();
//   int count = 0;
  for (; itlists != lists.end(); itlists ++, itposlists ++){
    std::list<PositionParticle * >::iterator itpart = itlists->begin();
    std::vector<std::vector<value_type > >::iterator itpos = itposlists->begin();
//     std::cout << count ++ << '\t' 
// 	      << itlists->size() << '\t' 
// 	      << itposlists->size() << std::endl;
    for (; itpart != itlists->end(); ++itpart, ++itpos){
      std::vector<value_type > diff ((*itpart)->r());
      VectorOperation::add (diff, -1, *itpos);
      pbox->shortestImage (diff);
      if (VectorOperation::dot (diff, diff) > value){
	return true;
      }
    }
  }  
  return false;
}

void NeighborListPairFinder::checkList ()
{
  if (checkMove()){
//     std::cout << "rebuild!!\n";
    rebuild();
  }
}

void NeighborListPairFinder::unregisterParticle (PositionParticle & p)
{
  for (std::vector<std::list<PositionParticle * > >::iterator itLists = lists.begin();
       itLists != lists.end(); itLists ++){
    for (std::list<PositionParticle * >::iterator itParts = itLists->begin();
	 itParts != itLists->end(); itParts ++){
      if ((*itParts)->id() == p.id()){
	itLists->erase (itParts);
	// std::cout << "deleted\n";
	return;
      }
    }
  }
}

void NeighborListPairFinder::unregisterParticle (std::vector<PositionParticle * > & p)
{  
  for (std::vector<PositionParticle * >::iterator itp = p.begin();
       itp != p.end(); ){
    for (std::vector<std::list<PositionParticle * > >::iterator itLists = lists.begin();
	 itLists != lists.end(); itLists ++){
      for (std::list<PositionParticle * >::iterator itParts = itLists->begin();
	   itParts != itLists->end(); itParts ++){
	if ((*itParts)->id() == (*itp)->id()){
	  itLists->erase (itParts);
	  std::cout << "deleted " << (*itp)->id() << std::endl;
	  goto loop; 
	}
      }      
    }
    loop:
    itp ++;
  }
}

void makeSwitchVector (std::vector<std::vector<int > > & vecs)
{
  vecs.resize(14);
  std::vector<int > tmp(3);
  int i = 0;
  tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
  vecs[i++] = tmp;
  tmp[0] = 1; tmp[1] = 0; tmp[2] = 0;
  vecs[i++] = tmp;
  tmp[0] = 1; tmp[1] = 1; tmp[2] = 0;
  vecs[i++] = tmp;
  tmp[0] = 0; tmp[1] = 1; tmp[2] = 0;
  vecs[i++] = tmp;
  tmp[0] = -1; tmp[1] = 1; tmp[2] = 0;
  vecs[i++] = tmp;
  tmp[0] = 0; tmp[1] = 0; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = 1; tmp[1] = 0; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = 1; tmp[1] = 1; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = 0; tmp[1] = 1; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = -1; tmp[1] = 1; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = -1; tmp[1] = 0; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = -1; tmp[1] = -1; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = 0; tmp[1] = -1; tmp[2] = 1;
  vecs[i++] = tmp;
  tmp[0] = 1; tmp[1] = -1; tmp[2] = 1;
  vecs[i++] = tmp;
}


// inline bool AllPairPairFinder::getPair (pPositionParticle & p0, 
// 					pPositionParticle & p1)
//   if (iter1 != list.end()){
//     p0 = (*iter0);
//     p1 = (*iter1);
//     iter1 ++;
//     return true;
//   }
//   else {
//     iter0 ++;
//     preiter0 ++;
//     if (preiter0 == list.end()){
//       return false;
//     }
//     else {
//       iter1 = preiter0;
//       p0 = (*iter0);
//       p1 = (*iter1);
//       iter1 ++;
//       return true;
//     }
//   }

      
