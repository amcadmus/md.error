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

void loadFrame (const char * filename,
		std::vector<std::vector<double > > &coord,
		std::vector<std::vector<double > > &force)
{
  FILE * fp = fopen (filename, "r");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", filename);
    exit (1);
  }  
  unsigned i = 0;
  
  while (fscanf (fp, "%lf %lf %lf %lf %lf %lf",
		 &coord[i][0],
		 &coord[i][1],
		 &coord[i][2],
		 &force[i][0],
		 &force[i][1],
		 &force[i][2]) == 6){
    i ++;
  }

  fclose (fp);
}

// void loadFrame (const char * filename,
// 		CoordType * coord,
// 		ScalorType * forcx,
// 		ScalorType * forcy,
// 		ScalorType * forcz)
// {
//   FILE * fp = fopen (filename, "r");
//   if (fp == NULL){
//     fprintf (stderr, "cannot open file %s", filename);
//     exit (1);
//   }  
//   unsigned i = 0;
  
//   while (fscanf (fp, "%f %f %f %f %f %f",
// 		 &coord[i].x,
// 		 &coord[i].y,
// 		 &coord[i].z,
// 		 &forcx[i],
// 		 &forcy[i],
// 		 &forcz[i]) != 6){
//     i ++;
//   }

//   fclose (fp);
// }



int main(int argc, char * argv[])
{
  char * filename;
  char dirname [1024];

  if (argc != 6){
    printf ("Usage:\n%s conf.gro dirName device refh rcut2\n", argv[0]);
    return 1;
  }
  filename = argv[1];
  strncpy (dirname, argv[2], 1024);
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

  std::vector<std::vector<double > > coord, force;
  coord.resize (sys.hdata.numAtom, std::vector<double > (3, 0.));
  force.resize (sys.hdata.numAtom, std::vector<double > (3, 0.));
  // CoordType * stored_coord;
  // ScalorType * stored_forcex;
  // ScalorType * stored_forcey;
  // ScalorType * stored_forcez;
  // stored_coord = (CoordType *) malloc (sizeof(CoordType) * sys.hdata.numAtom);
  // stored_forcex = (ScalorType *) malloc (sizeof(ScalorType) * sys.hdata.numAtom);
  // stored_forcey = (ScalorType *) malloc (sizeof(ScalorType) * sys.hdata.numAtom);
  // stored_forcez = (ScalorType *) malloc (sizeof(ScalorType) * sys.hdata.numAtom);

  std::vector<double > boxsize (3);
  boxsize[0] = sys.box.size.x;
  boxsize[1] = sys.box.size.y;
  boxsize[2] = sys.box.size.z;

  ErrorProfile_PiecewiseConst ep (boxsize, refh);

  char timename[1024];
  sprintf (timename, "%s/times.out", dirname);
  FILE * fptime = fopen (timename, "r");
  if (fptime == NULL){
    fprintf (stderr, "cannot open file %s\n", timename);
    exit (1);
  }
  ScalorType nowtime;
  while (fscanf (fptime, "%f", &nowtime) == 1){
    printf ("loaded frame at time %f ps       \r", nowtime);
    fflush (stdout);
    char nowfile [1024];
    sprintf (nowfile, "%s/posiForc_t%05.3f.out", dirname);
    // loadFrame (nowfile, stored_coord, stored_forcex, stored_forcey, stored_forcez);
    loadFrame (nowfile, coord, force);
    
    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      sys.hdata.coord[i].x = coord[i][0];
      sys.hdata.coord[i].y = coord[i][1];
      sys.hdata.coord[i].z = coord[i][2];
      // sys.hdata.coord[i].x = stored_coord[i].x;
      // sys.hdata.coord[i].y = stored_coord[i].y;
      // sys.hdata.coord[i].z = stored_coord[i].z;
    }
    cpyHostMDDataToDevice (&sys.hdata, &sys.ddata);
    inter.clearInteraction (sys);
    inter.applyNonBondedInteraction (sys, rcut2);
    cpyDeviceMDDataToHost (&sys.ddata, &sys.hdata);

    for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
      // coord[i][0] = sys.hdata.coord[i].x;
      // coord[i][1] = sys.hdata.coord[i].y;
      // coord[i][2] = sys.hdata.coord[i].z;
      // stored_forcex[i] -= sys.hdata.forcx[i];
      // stored_forcey[i] -= sys.hdata.forcy[i];
      // stored_forcez[i] -= sys.hdata.forcz[i];
      force[i][0] -= sys.hdata.forcx[i];
      force[i][1] -= sys.hdata.forcy[i];
      force[i][2] -= sys.hdata.forcz[i];
      double tmp = sqrt (force[i][0]*force[i][0] +
			 force[i][1]*force[i][1] +
			 force[i][2]*force[i][2] );
      if (tmp > 0.01){
	printf ("time: %f, tmp is %f, i is %d, coordx is %f\n", nowtime, tmp, i, coord[i][0]);
      }
    }
    ep.deposit (coord, force);
  }
  ep.calculate();
  ep.print_x (("real.x.out"));
  ep.print_xy (("real.xy.out"));
  
  // free (stored_coord);
  // free (stored_forcex);
  // free (stored_forcey);
  // free (stored_forcez);
  
  return 0;
}

  
