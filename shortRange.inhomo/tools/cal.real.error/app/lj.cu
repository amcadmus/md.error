#include <stdio.h>
#include "MDSystem_interface.h"
#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "RandomGenerator.h"
#include "Auxiliary.h"
#include "NeighborList_interface.h"
#include"Statistic.h"
#include "Integrator_interface.h"
#include "InteractionEngine_interface.h"
#include "tmp.h"
#include "Reshuffle_interface.h"
#include "Displacement_interface.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
  ScalorType rcut1 = 5.f;
  char * filename;
  
  if (argc != 4){
    printf ("Usage:\n%s conf.gro rcut1 device\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    rcut1 = atof(argv[2]);
    filename = argv[1];
  }
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");

  MDSystem sys;
  sys.initConfig(filename);

  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  ScalorType rcut2 = sys.box.size.z / 2 - 1.f;
  ljparam.reinit (1.f, 1.f, 0.f, rcut1, rcut2);
  
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = rcut2;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, 10.f);
  sys.normalizeDeviceData ();
  clist.rebuild (sys, NULL);
  nlist.rebuild (sys, clist, NULL);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  
  try{ 
    inter.clearInteraction (sys);
    inter.applyNonBondedInteraction (sys, nlist, NULL, NULL);

    sys.updateHostFromDevice (NULL);
    FILE *fp = fopen ("force.out", "w");
    fprintf (fp, "%d\n%f %f %f\n",
	     sys.ddata.numAtom,
	     sys.box.size.x, sys.box.size.y, sys.box.size.z);
    for (unsigned i = 0; i < sys.ddata.numAtom; ++i){
      fprintf (fp, "%e %e %e  %e %e %e\n",
	       sys.hdata.coord[i].x, 
	       sys.hdata.coord[i].y, 
	       sys.hdata.coord[i].z,
	       sys.hdata.forcx[i],
	       sys.hdata.forcy[i],
	       sys.hdata.forcz[i]);
    }
    fclose (fp);
  }
  catch (MDException &e){
    fprintf (stderr, "%s\n", e.what());
    return 1;
  }
  
  
  return 0;
}

  
