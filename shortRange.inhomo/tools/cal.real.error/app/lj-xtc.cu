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

#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
  ScalorType rcut1 = 5.f;
  char * filename;
  char * xtcname;

  if (argc != 6){
    printf ("Usage:\n%s conf.gro traj.xtc rcut1 device refh\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    rcut1 = atof(argv[3]);
    filename = argv[1];
    xtcname  = argv[2];
  }
  printf ("# setting device to %d\n", atoi(argv[4]));
  cudaSetDevice (atoi(argv[4]));
  checkCUDAError ("set device");
  double refh = atof (argv[5]);

  MDSystem sys;
  sys.initConfig(filename);

  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  // ScalorType rcut2 = sys.box.size.z / 2 - 1.f;
  ScalorType rcut2 = rcut1 * 3;
  if (rcut2 > sys.box.size.z / 2.f - 1.f) rcut2 = sys.box.size.z / 2 - 1.f;
  printf ("# rcut1 is %f\n", rcut1);
  printf ("# rcut2 is %f\n", rcut2);
  ljparam.reinit (1.f, 1.f, 0.f, rcut1, rcut2);
  
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  // ScalorType maxrcut = sysNbInter.maxRcut();
  // ScalorType rlist = rcut2;
  // CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  // NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, 2.f);
  // sys.normalizeDeviceData ();
  // clist.rebuild (sys, NULL);
  // nlist.rebuild (sys, clist, NULL);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  

  inter.clearInteraction (sys);
  inter.applyNonBondedInteraction (sys, rcut2);
    
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

  
  
  int step;
  int natoms= 0;
  float time, prec;
  matrix gbox;  
  rvec * xx;
  xx = (rvec *) malloc (sizeof(rvec) * sys.hdata.numAtom);
  XDRFILE * fpxtc = xdrfile_open (xtcname, "r");
  if (fpxtc == NULL){
    fprintf (stderr, "cannot open file %s\n", xtcname);
    return 1;;
  }
  std::vector<double > boxsize (3);
  boxsize[0] = sys.box.size.x;
  boxsize[1] = sys.box.size.y;
  boxsize[2] = sys.box.size.z;
  ErrorProfile_PiecewiseConst ep (boxsize, refh);
  std::vector<std::vector<double > > coord, force;
  coord.resize (sys.hdata.numAtom, std::vector<double > (3, 0.));
  force.resize (sys.hdata.numAtom, std::vector<double > (3, 0.));
  while (read_xtc (fpxtc, natoms, &step, &time, gbox, xx, &prec) == 0){
    printf ("loaded frame at time %f ps       \r", time);
    fflush (stdout);
    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      sys.hdata.coord[i].x = xx[i][0];
      sys.hdata.coord[i].y = xx[i][1];
      sys.hdata.coord[i].z = xx[i][2];
    }
    cpyHostMDDataToDevice (&sys.hdata, &sys.ddata);
    inter.clearInteraction (sys);
    inter.applyNonBondedInteraction (sys, rcut2);
    cpyDeviceMDDataToHost (&sys.ddata, &sys.hdata);

    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      coord[i][0] = sys.hdata.coord[i].x;
      coord[i][1] = sys.hdata.coord[i].y;
      coord[i][2] = sys.hdata.coord[i].z;
      force[i][0] = sys.hdata.forcx[i];
      force[i][1] = sys.hdata.forcy[i];
      force[i][2] = sys.hdata.forcz[i];
    }
    ep.deposit (coord, force);
  }
  ep.calculate();
  ep.print_x (("real.x.out"));
  ep.print_xy (("real.xy.out"));
  xdrfile_close (fpxtc);
  free (xx);
  
  return 0;
}

  
