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
#include "AssignForceCorr.h"

#define NThreadsPerBlockCell	128
#define NThreadsPerBlockAtom	96

int main(int argc, char * argv[])
{
  char * filename;
  char xtcname[] = "traj.xtc";
  char afcname[] = "fc.ctj";

  if (argc != 7){
    printf ("Usage:\n%s conf.gro device refh rcut1 rcut2 start_t\n", argv[0]);
    return 1;
  }
  filename = argv[1];
  printf ("# setting device to %d\n", atoi(argv[2]));
  cudaSetDevice (atoi(argv[2]));
  checkCUDAError ("set device");
  double refh = atof (argv[3]);
  ScalorType start_t = atof (argv[6]);
  
  MDSystem sys;
  sys.initConfig(filename);

  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  ScalorType rcut1 = atof(argv[4]);
  ScalorType rcut2 = atof(argv[5]);
  printf ("# rcut2 is %f\n", rcut2);
  int nimage = (rcut2 - 0.00001) / sys.box.size.y;
  nimage ++;
  printf ("#@ nimage is %d\n", nimage);

  ljparam.reinit (1.f, 1.f, 0.f, rcut1, rcut2);
  
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);

  AssignForceCorr afc;
  afc.reinit (sys, NThreadsPerBlockAtom);
  afc.init_read (afcname);
  
  int step;
  int natoms= 0;
  float time, afctime, prec;
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
    afc.read (afctime);
    if (fabs (time - afctime) > 1e-4) {
      printf ("inconsistent trajactories\n");
      exit (1);
    }
    if (time < start_t - 1e-4) continue;
    printf ("loaded frame at time %f ps       \r", time);
    fflush (stdout);
    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      sys.hdata.coord[i].x = xx[i][0];
      sys.hdata.coord[i].y = xx[i][1];
      sys.hdata.coord[i].z = xx[i][2];
    }
    cpyHostMDDataToDevice (&sys.hdata, &sys.ddata);
    inter.clearInteraction (sys);
    inter.applyNonBondedInteraction (sys, rcut2, NULL);
    cpyDeviceMDDataToHost (&sys.ddata, &sys.hdata);
    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      coord[i][0] = sys.hdata.coord[i].x;
      coord[i][1] = sys.hdata.coord[i].y;
      coord[i][2] = sys.hdata.coord[i].z;
      ScalorType cfx, cfy, cfz;
      afc.getForceCorr (coord[i][0], coord[i][1], coord[i][2],
			cfx, cfy, cfz);
      force[i][0] = sys.hdata.forcx[i] - cfx;
      force[i][1] = sys.hdata.forcy[i] - cfy;
      force[i][2] = sys.hdata.forcz[i] - cfz;
    }
    ep.deposit (coord, force);
  }
  
  ep.calculate();
  // ep.print_x (("real.x.out"));
  ep.print_x_avg (("a.real.x.out"));
  // ep.print_xy (("real.xy.out"));
  
  return 0;
}

  
