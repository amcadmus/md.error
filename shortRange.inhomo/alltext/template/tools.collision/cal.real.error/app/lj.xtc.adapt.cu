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
  char xtcname[] = "traj.xtc";
  float start_t;
  
  if (argc != 7){
    printf ("Usage:\n%s conf.gro rcut.rtj device refh rcut2 start_t\n", argv[0]);
    return 1;
  }
  filename = argv[1];
  strncpy (rcutsavename, argv[2], 1024);
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");
  double refh = atof (argv[4]);
  start_t = atof (argv[6]);

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

  AssignRCut assign;
  assign.reinit (sys, NThreadsPerBlockAtom);
  assign.init_read (rcutsavename);

  int step;
  int natoms= 0;
  float time, timercut, prec;
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
  read_xtc (fpxtc, natoms, &step, &time, gbox, xx, &prec);
  while (read_xtc (fpxtc, natoms, &step, &time, gbox, xx, &prec) == 0){
    assign.read (timercut);
    if (fabs (time - timercut) > 1e-4) {
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
    assign.assign (sys);
    inter.clearInteraction (sys);
    inter.applyNonBondedInteraction (sys, NULL, rcut2);
    cpyDeviceMDDataToHost (&sys.ddata, &sys.hdata);
    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      coord[i][0] = sys.hdata.coord[i].x;
      coord[i][1] = sys.hdata.coord[i].y;
      coord[i][2] = sys.hdata.coord[i].z;
      force[i][0] = sys.hdata.forcx[i];
      force[i][1] = sys.hdata.forcy[i];
      force[i][2] = sys.hdata.forcz[i];
    }
    // {
    //   int i = 179;
    //   double fsumx, fsumy, fsumz;
    //   fsumx = fsumy = fsumz = 0.;
    //   for (int jj = 0; jj < sys.hdata.numAtom; ++jj){
    // 	double diffx = coord[i][0] - coord[jj][0];
    // 	if      (diffx < -0.5 * sys.box.size.x) diffx += sys.box.size.x;
    // 	else if (diffx >= 0.5 * sys.box.size.x) diffx -= sys.box.size.x;
    // 	for (int idy = -2; idy <= 2; ++idy){
    // 	  for (int idz = -2; idz <= 2; ++idz){
    // 	    double diffy = coord[i][1] - coord[jj][1] - idy * sys.box.size.y;
    // 	    double diffz = coord[i][2] - coord[jj][2] - idz * sys.box.size.z;
    // 	    if (idy == 0 && idz == 0 && jj == i) continue;
    // 	    ScalorType dr2 = diffx*diffx+diffy*diffy+diffz*diffz;
    // 	    if ((dr2 <= 30*30) && (dr2 >= 5*5)){
    // 	      double boolScalor = 4.f;
    // 	      double ri2 = 1.f/dr2;
    // 	      double sri2 = ri2;
    // 	      ScalorType sri6 = sri2*sri2*sri2;
    // 	      boolScalor *= (6.f * sri6 - 12.f * sri6*sri6) * ri2;
    // 	      fsumx += diffx * boolScalor;
    // 	      fsumy += diffy * boolScalor;
    // 	      fsumz += diffz * boolScalor;
    // 	    }
    // 	  }
    // 	}
    //   }
    //   printf ("%f\n", fsumx);
    // }
    ep.deposit (coord, force);
    ep.calculate();
    char timefile[1024];
    int zhengshu = int(time);
    int xiaoshu  = int((time - zhengshu) * 100.);
    sprintf (timefile, "error.t%04d.%02d.out", zhengshu, xiaoshu);
    ep.print_xy (timefile);
    ep.clear();
  }
  
  // ep.calculate();
  // ep.print_x (("real.x.out"));
  // ep.print_x_avg (("a.real.x.out"));
  // ep.print_xy (("real.xy.out"));  
  
  return 0;
}

  
