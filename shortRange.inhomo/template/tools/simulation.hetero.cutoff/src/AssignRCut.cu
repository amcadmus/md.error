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
reinit (const MDSystem & sys,
	const AdaptRCut & arc,
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
  
  nx = arc.getNx();
  ny = arc.getNy();
  nz = arc.getNz();
  nele = nx * ny * nz;
  box = sys.box;
  
  hrcut = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &drcut, sizeof(ScalorType ) * nele);
  checkCUDAError ("AssignRCut::reinit malloc drcut");
  malloced = true;
}

void AssignRCut::
getRCut (const AdaptRCut & arc)
{
  for (int i = 0; i < nele; ++i){
    hrcut[i] = arc.getRCut()[i];
  }
  cudaMemcpy (drcut, hrcut, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  checkCUDAError ("AssignRCut::getRCut copy");
}

void AssignRCut::
uniform (const double & rc)
{
  for (int i = 0; i < nele; ++i){
    hrcut[i] = rc;
  }
  cudaMemcpy (drcut, hrcut, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  checkCUDAError ("AssignRCut::getRCut copy");
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

      
void AssignRCut::    
print_x (const char * file) const 
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    // double sum = 0.;
    // for (int j = 0; j < ny; ++j){
    //   for (int k = 0; k < nz; ++k){
    // 	sum += profile[index3to1(i, j, k)];
    //   }
    // }
    fprintf (fp, "%f %e\n",
	     (i + 0.5) * box.size.x / double(nx),
	     hrcut [index3to1(i, 0, 0)]
	);
  }
  fclose (fp);
}

void AssignRCut::
init_write (const char * file) const
{
  fp_write = fopen (file, "w");
  if (fp_write == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit(1);
  }
  double tmpbox[3];
  tmpbox[0] = box.size.x;
  tmpbox[1] = box.size.y;
  tmpbox[2] = box.size.z;
  int tmpnn[3];
  tmpnn[0] = nx;
  tmpnn[1] = ny;
  tmpnn[2] = nz;
  
  fwrite (tmpbox, sizeof(double), 3, fp_write);
  fwrite (tmpnn,  sizeof(int),    3, fp_write);
}

void AssignRCut::
end_write () const
{
  fclose (fp_write);
}

void AssignRCut::
write (const ScalorType & time) const
{
  ScalorType tmptime = time;
  fwrite (&tmptime, sizeof(ScalorType), 1,    fp_write);
  fwrite (hrcut,    sizeof(ScalorType), nele, fp_write);
}



