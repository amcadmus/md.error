#include "AssignRCut.h"
#include "BoxGeometry.h"

AssignRCut::
AssignRCut ()
    : malloced (false)
{
}

AssignRCut::
~AssignRCut ()
{
  freeAll();
}


void AssignRCut::
reinit (const char * filename,
	const MDSystem & sys,
    	const IndexType & NThread)
{
  freeAll ();
  
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  FILE * fp = fopen (filename, "r");  
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", filename);
    exit(1);    
  }
  fscanf (fp, "%d %d %d\n", &nx, &ny, &nz);
  nele = nx * ny * nz;
  box = sys.box;
  // hx = sys.box.size.x / nx;
  // hy = sys.box.size.y / ny;
  // hz = sys.box.size.z / nz;
  
  hrcut = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &drcut, sizeof(ScalorType ) * nele);
  checkCUDAError ("AssignRCut::reinit malloc drcut");
  malloced = true;

  maxRCut = 0.;
  for (int i = 0; i < nele; ++i){
    int c;
    c = fscanf (fp, "%f", &hrcut[i]);
    if (c != 1){
      printf ("c is not 1\n");
      exit (1);
    }
    if (hrcut[i] > maxRCut) maxRCut = hrcut[i];
  }
  cudaMemcpy (drcut, hrcut, sizeof(ScalorType ) * nele, cudaMemcpyHostToDevice);
  checkCUDAError ("AssignRCut::reinit cpy rcut");
  
  fclose (fp);
}

void AssignRCut::
freeAll ()
{
  if (malloced) {
    cudaFree (drcut);
    free (hrcut);
    malloced = false;
  }
}

using namespace RectangularBoxGeometry;

static void __global__
assignRCutToSystem (const ScalorType * rcutLattice,
		    const int nx,
		    const int ny,
		    const int nz,
		    const RectangularBox box,
		    const CoordType * coord,
		    const int numAtom,
		    ScalorType * rcut)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom){
    CoordType mycoord = coord[ii];
    if (mycoord.x < 0) mycoord.x += box.size.x;
    else if (mycoord.x >= box.size.x) mycoord.x -= box.size.x;
    if (mycoord.y < 0) mycoord.y += box.size.y;
    else if (mycoord.y >= box.size.y) mycoord.y -= box.size.y;
    if (mycoord.z < 0) mycoord.z += box.size.z;
    else if (mycoord.z >= box.size.z) mycoord.z -= box.size.z;
    int ix = (mycoord.x * nx) / box.size.x;
    int iy = (mycoord.y * ny) / box.size.y;
    int iz = (mycoord.z * nz) / box.size.z;
    int idx = iz + nz * (iy + ny * ix);
    rcut[ii] = rcutLattice[idx];
  }
}

void AssignRCut::
assign (MDSystem & sys)
{
  assignRCutToSystem
      <<<atomGridDim, myBlockDim>>> (
	  drcut,
	  nx, ny, nz,
	  sys.box,
	  sys.ddata.coord,
	  sys.ddata.numAtom,
	  sys.ddata.rcut);
}

      
