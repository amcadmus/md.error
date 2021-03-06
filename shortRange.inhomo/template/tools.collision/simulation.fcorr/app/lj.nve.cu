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
#include "AssignForceCorr.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fftw3.h>

// #define NThreadsPerBlockCell	32
// #define NThreadsPerBlockAtom	4

#define NThreadsPerBlockCell	96
#define NThreadsPerBlockAtom	96

#include "DensityProfile.h"

int main(int argc, char * argv[])
{
  IndexType nstep = 100000;
  IndexType confFeq = 2000;
  IndexType thermoFeq = 100;
  ScalorType dt = 0.002;
  ScalorType rcut = 2.5;
  ScalorType nlistExten = 0.49;
  ScalorType nlistSizeFactor = 100.0;
  char * filename;

  double refh = 1.0;
  IndexType densityProfileSamplingFeq = 10;
  IndexType assignForceCorrFeq = 10;
  IndexType updateForceCorrFeq = 10;
  char afcname[] = "fc.ctj";

  if (argc != 4){
    printf ("Usage:\n%s conf.gro nstep device\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    nstep = atoi(argv[2]);
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
  ljparam.reinit (1.f, 1.f, 0.f, rcut);
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  DensityProfile_PiecewiseConst dp;
  printf ("# init DensityProfile_PiecewiseConst\n");
  dp.reinit (sys.box.size.x, sys.box.size.y, sys.box.size.z, refh);
  ForceCorr fc;
  printf ("# init AdaptRCut\n");
  fc.reinit (rcut, dp);
  AssignForceCorr afc;
  printf ("# init AssignRCut\n");
  afc.reinit (sys, fc, NThreadsPerBlockAtom);
  afc.zero ();
  afc.print_x ("fc.x.out");
  afc.assign (sys);
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  ScalorType energyCorr = sysNbInter.energyCorrection ();
  ScalorType pressureCorr = sysNbInter.pressureCorrection ();
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = maxrcut + nlistExten;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, nlistSizeFactor);
  sys.normalizeDeviceData ();
  clist.rebuild (sys, NULL);
  nlist.rebuild (sys, clist, NULL);
  Displacement_max disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  
  MDStatistic st(sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  
  MDTimer timer;
  unsigned i;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);

  Reshuffle resh (sys);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
  }
  
  printf ("# prepare ok, start to run\n");
  sys.recoverDeviceData (&timer);
  sys.updateHostFromRecovered (&timer);
  sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6       7\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  Npairs\n");

  try{
    sys.initWriteXtc ("traj.xtc");
    sys.initWriteTrr ("traj.trr");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    sys.writeHostDataTrr (0, 0*dt, &timer);
    afc.init_write (afcname);
    for (i = 0; i < nstep; ++i){
      if (i%1 == 0){
	tfremover.remove (sys, &timer);
      }
      
      inte_vv.step1 (sys, dt, &timer);

      st.clearDevice();
      inter.clearInteraction (sys);
      ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
      if (maxdr > nlistExten * 0.5){
	// printf ("# Rebuild at step %09i ... ", i+1);
	// fflush(stdout);
	// rebuild
	sys.normalizeDeviceData (&timer);
	disp.recordCoord (sys);
	clist.rebuild (sys, &timer);
	nlist.rebuild (sys, clist, &timer);
	// printf ("done\n");
	// fflush(stdout);
      }
      if ((i+1) % thermoFeq == 0){	
	inter.applyNonBondedInteraction (sys, nlist, st, NULL, &timer);
      }
      else {
	inter.applyNonBondedInteraction (sys, nlist, NULL, &timer);
      }      

      if ((i) % assignForceCorrFeq == 0){
        afc.assign (sys);
      }
      
      if ((i+1) % thermoFeq == 0){	
	inte_vv.step2 (sys, dt, st, &timer);
      }
      else {
	inte_vv.step2 (sys, dt, &timer);
      }      

      if ((i+1) % thermoFeq == 0){
	timer.tic (mdTimeDataIO);
	st.updateHost ();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.5e\n",
		(i+1),  
		(i+1) * dt, 
		st.nonBondedEnergy(),
		st.kineticEnergy(),
		st.kineticEnergy() * 2. / 3. / (double (sys.hdata.numAtom) - 3.),
		st.nonBondedEnergy() +
		st.kineticEnergy(),
		double (nlist.calSumNeighbor ())
	    );
	fflush(stdout);
	timer.toc (mdTimeDataIO);
      }

      if ((i+1) % densityProfileSamplingFeq == 0) {
	timer.tic (mdTimeDensityProfile);
	sys.updateHostFromDevice (NULL);
	dp.deposite (sys.hdata.coord, sys.hdata.numAtom);
	timer.toc (mdTimeDensityProfile);
      }

      if ((i+1) % updateForceCorrFeq == 0) {
	// printf ("# update rcut\n");
	timer.tic (mdTimeDensityProfile);
	dp.calculate ();
	dp.print_x ("density.x.out");
	timer.toc (mdTimeDensityProfile);
	timer.tic (mdTimeAdaptRCut);
	fc.calError (dp);
	// fc.print_x ("fc.x.out");
	afc.getForceCorr (fc);
	afc.print_x ("afc.x.out");
	timer.toc (mdTimeAdaptRCut);
	if (i != nstep - 1) dp.clearData ();
      }
      
      if ((i+1) % confFeq == 0){
      	sys.recoverDeviceData (&timer);
      	sys.updateHostFromRecovered (&timer);
      	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      	sys.writeHostDataTrr (i+1, (i+1)*dt, &timer);
	afc.write ((i+1)*dt);
      }
      
      if ((i+1) % 100 == 0){
	if (resh.calIndexTable (clist, &timer)){
	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
	  clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
	}
      }
    }
    sys.endWriteXtc();
    sys.endWriteTrr();
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataGro ("confout.gro", nstep, nstep*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
  }
  catch (MDExcptCuda & e){
    // resh.recoverMDDataToHost (sys, &timer);
    // sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
    return 1;
  }
  catch (MDException &e){
    fprintf (stderr, "%s\n", e.what());
    return 1;
  }

  // dp.end_write();
  dp.save ("density.save");

  afc.end_write();
  
  return 0;
}

  
