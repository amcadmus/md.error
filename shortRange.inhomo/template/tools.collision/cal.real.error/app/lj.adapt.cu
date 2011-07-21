#define CPLUSPLUS

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

#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
// #include "GroFileManager.h"
#include "ErrorProfile.h"
#include "AdaptRCut.h"
#include "AssignRCut.h"

#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
  char * filename;
  char rcutsavename [1024];

  if (argc != 6){
    printf ("Usage:\n%s conf.gro rcut.save device refh rcut2\n", argv[0]);
    return 1;
  }
  filename = argv[1];
  strncpy (rcutsavename, argv[2], 1024);
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");
  double refh = atof (argv[4]);

  MDSystem sys;
  sys.initConfig(filename);

  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  ScalorType rcut2 = atof(argv[5]);
  printf ("# rcut2 is %f\n", rcut2);
  int nimage = (rcut2 - 0.00001) / sys.box.size.y;
  nimage ++;
  printf ("#@ nimage is %d\n", nimage);

  ljparam.reinit (1.f, 1.f, 0.f, 0.f, rcut2);
  
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);

  inter.clearInteraction (sys);
  inter.applyNonBondedInteraction (sys, rcut2);    

  std::vector<double > boxsize (3);
  boxsize[0] = sys.box.size.x;
  boxsize[1] = sys.box.size.y;
  boxsize[2] = sys.box.size.z;

  ErrorProfile_PiecewiseConst ep (boxsize, refh);
  AdaptRCut arc;
  arc.load_rc (std::string(rcutsavename));
  AssignRCut assign;
  assign.reinit (sys, arc, NThreadsPerBlockAtom);
  assign.getRCut (arc);
  assign.assign (sys);

  inter.clearInteraction (sys);
  inter.applyNonBondedInteraction (sys, NULL, rcut2);
  cpyDeviceMDDataToHost (&sys.ddata, &sys.hdata);

  std::vector<std::vector<double > > coord, force;

  for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
    std::vector<double > tmp(3);
    tmp[0] = sys.hdata.coord[i].x;
    tmp[1] = sys.hdata.coord[i].y;
    tmp[2] = sys.hdata.coord[i].z;
    coord.push_back (tmp);
    tmp[0] = sys.hdata.forcx[i];
    tmp[1] = sys.hdata.forcy[i];
    tmp[2] = sys.hdata.forcz[i];
    force.push_back (tmp);
  }

  ep.deposit (coord, force);
  ep.calculate();
  // ep.print_x (("real.x.out"));
  ep.print_x_avg (("a.real.x.out"));
  // ep.print_xy (("real.xy.out"));
  
  return 0;
}

  
