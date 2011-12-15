#include "AssignForceCorr.h"
#include "BoxGeometry.h"

AssignForceCorr::
AssignForceCorr ()
    : malloced (false)
{
}

AssignForceCorr::
~AssignForceCorr ()
{
  freeAll();
}

void AssignForceCorr::
reinit (const MDSystem & sys,
	const ForceCorr & arc,
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
  
  hfcx = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  hfcy = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  hfcz = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  for (int i = 0; i < nele; ++i){
    hfcx[i] = hfcy[i] = hfcz[i] = 0.;
  }
  cudaMalloc ((void **) &dfcx, sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &dfcy, sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &dfcz, sizeof(ScalorType ) * nele);
  checkCUDAError ("AssignForceCorr::reinit malloc drcut");
  malloced = true;
}

__global__ void
copyForceCorr (const cufftComplex * error2rx,
	       const cufftComplex * error2ry,
	       const cufftComplex * error2rz,
	       const int nele,
	       ScalorType * dfcx,
	       ScalorType * dfcy,
	       ScalorType * dfcz)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;

  if (ii < nele){
    dfcx[ii] = error2rx[ii].x;
    dfcy[ii] = error2ry[ii].x;
    dfcz[ii] = error2rz[ii].x;
  }
}

__global__ void
zero_d (const int nele,
	ScalorType * dfcx,
	ScalorType * dfcy,
	ScalorType * dfcz)
{
  unsigned bid = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned tid = threadIdx.x;
  unsigned ii = tid + bid * blockDim.x;

  if (ii < nele){
    dfcx[ii] = 0.f;
    dfcy[ii] = 0.f;
    dfcz[ii] = 0.f;
  }
}


void AssignForceCorr::
getForceCorr (const ForceCorr & arc)
{
  unsigned blockSize = 128;
  unsigned nblock = unsigned(nele) / blockSize + 1;
  copyForceCorr <<<nblock, blockSize>>>
      (arc.d_error2rx,
       arc.d_error2ry,
       arc.d_error2rz,
       nele,
       dfcx,
       dfcy,
       dfcz);
  // for (int i = 0; i < nele; ++i){
  //   hfcx[i] = arc.getForceCorrX ()[i][0];
  //   hfcy[i] = arc.getForceCorrY ()[i][0];
  //   hfcz[i] = arc.getForceCorrZ ()[i][0];    
  // }
  // cudaMemcpy (dfcx, hfcx, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  // cudaMemcpy (dfcy, hfcy, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  // cudaMemcpy (dfcz, hfcz, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  checkCUDAError ("AssignForceCorr::getRCut copy");
}

void AssignForceCorr::
zero ()
{
  unsigned blockSize = 128;
  unsigned nblock = unsigned(nele) / blockSize + 1;
  zero_d <<<nblock, blockSize>>> (nele, dfcx, dfcy, dfcz);
  // for (int i = 0; i < nele; ++i){
  //   hfcx[i] = 0.;
  //   hfcy[i] = 0.;
  //   hfcz[i] = 0.;
  // }
  // cudaMemcpy (dfcx, hfcx, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  // cudaMemcpy (dfcy, hfcy, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  // cudaMemcpy (dfcz, hfcz, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  checkCUDAError ("AssignForceCorr::getRCut copy");
}



void AssignForceCorr::
freeAll ()
{
  if (malloced) {
    cudaFree (dfcx);
    cudaFree (dfcy);
    cudaFree (dfcz);
    free (hfcx);
    free (hfcy);
    free (hfcz);
    malloced = false;
  }
}

using namespace RectangularBoxGeometry;

static void __global__
assignForceToSystem (const ScalorType * dfcx,
		     const ScalorType * dfcy,
		     const ScalorType * dfcz,
		     const int nx,
		     const int ny,
		     const int nz,
		     const RectangularBox box,
		     const CoordType * coord,
		     const int numAtom,
		     ScalorType * fcx,
		     ScalorType * fcy,
		     ScalorType * fcz)
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
    fcx[ii] = dfcx[idx];
    fcy[ii] = dfcy[idx];
    fcz[ii] = dfcz[idx];
  }
}

void AssignForceCorr::
assign (MDSystem & sys)
{
  assignForceToSystem
      <<<atomGridDim, myBlockDim>>> (
	  dfcx, dfcy, dfcz,
	  nx, ny, nz,
	  sys.box,
	  sys.ddata.coord,
	  sys.ddata.numAtom,
	  sys.ddata.fcx,
	  sys.ddata.fcy,
	  sys.ddata.fcz
	  );
}

      
void AssignForceCorr::    
print_x (const char * file) const 
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }
  cudaMemcpy (hfcx, dfcx, sizeof(ScalorType)*nele, cudaMemcpyDeviceToHost);
  cudaMemcpy (hfcy, dfcy, sizeof(ScalorType)*nele, cudaMemcpyDeviceToHost);
  cudaMemcpy (hfcz, dfcz, sizeof(ScalorType)*nele, cudaMemcpyDeviceToHost);

  for (int i = 0; i < nx; ++i){
    // double sum = 0.;
    // for (int j = 0; j < ny; ++j){
    //   for (int k = 0; k < nz; ++k){
    // 	sum += profile[index3to1(i, j, k)];
    //   }
    // }
    fprintf (fp, "%f %e %e %e\n",
	     (i + 0.5) * box.size.x / double(nx),
	     hfcx [index3to1(i, 0, 0)],
	     hfcy [index3to1(i, 0, 0)],
	     hfcz [index3to1(i, 0, 0)]
	);
  }
  fclose (fp);
}

void AssignForceCorr::
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

void AssignForceCorr::
end_write () const
{
  fclose (fp_write);
}

void AssignForceCorr::
write (const ScalorType & time) const
{
  ScalorType tmptime = time;
  cudaMemcpy (hfcx, dfcx, sizeof(ScalorType)*nele, cudaMemcpyDeviceToHost);
  cudaMemcpy (hfcy, dfcy, sizeof(ScalorType)*nele, cudaMemcpyDeviceToHost);
  cudaMemcpy (hfcz, dfcz, sizeof(ScalorType)*nele, cudaMemcpyDeviceToHost);
  fwrite (&tmptime, sizeof(ScalorType), 1,    fp_write);
  fwrite (hfcx,     sizeof(ScalorType), nele, fp_write);
  fwrite (hfcy,     sizeof(ScalorType), nele, fp_write);
  fwrite (hfcz,     sizeof(ScalorType), nele, fp_write);
}



