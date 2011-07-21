#include "AssignForceCorr.h"
#include "BoxGeometry.h"
#include <iostream>

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
}


void AssignForceCorr::
zero ()
{
  for (int i = 0; i < nele; ++i){
    hfcx[i] = 0.;
    hfcy[i] = 0.;
    hfcz[i] = 0.;
  }
  cudaMemcpy (dfcx, hfcx, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  cudaMemcpy (dfcy, hfcy, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  cudaMemcpy (dfcz, hfcz, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
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


      
void AssignForceCorr::    
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
  fwrite (&tmptime, sizeof(ScalorType), 1,    fp_write);
  fwrite (hfcx,     sizeof(ScalorType), nele, fp_write);
  fwrite (hfcy,     sizeof(ScalorType), nele, fp_write);
  fwrite (hfcz,     sizeof(ScalorType), nele, fp_write);
}


void AssignForceCorr::
init_read (const char * file)
{
  fp_read = fopen (file, "r");
  if (fp_read == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit(1);
  }
  double tmpbox[3];
  int tmpnn[3];
  int status;
  
  status = fread (tmpbox, sizeof(double), 3, fp_read);
  if (status != 3){
    printf ("format error\n");
    exit(1);
  }
  status = fread (tmpnn,  sizeof(int),    3, fp_read);
  if (status != 3){
    printf ("format error\n");
    exit(1);
  }

  nx = tmpnn[0];
  ny = tmpnn[1];
  nz = tmpnn[2];
  nele = nx * ny * nz;
  setBoxSize (tmpbox[0], tmpbox[1], tmpbox[2], box);

  freeAll();
  hfcx = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  hfcy = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  hfcz = (ScalorType *) malloc (sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &dfcx, sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &dfcy, sizeof(ScalorType ) * nele);
  cudaMalloc ((void **) &dfcz, sizeof(ScalorType ) * nele);
  checkCUDAError ("AssignForceCorr::reinit malloc drcut");
  malloced = true;
}

void AssignForceCorr::
read (ScalorType & time) 
{
  int status;
  status = fread (&time, sizeof(ScalorType), 1,    fp_read);
  if (status != 1){
    printf ("format error\n");
    exit(1);
  }
  status = fread (hfcx, sizeof(ScalorType), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    exit(1);
  }
  status = fread (hfcy, sizeof(ScalorType), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    exit(1);
  }
  status = fread (hfcz, sizeof(ScalorType), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    exit(1);
  }
  cudaMemcpy (dfcx, hfcx, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  cudaMemcpy (dfcy, hfcy, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  cudaMemcpy (dfcz, hfcz, sizeof(ScalorType) * nele, cudaMemcpyHostToDevice);
  checkCUDAError ("AssignForceCorr::getRCut copy");
}

void AssignForceCorr::
end_read () const
{
  fclose (fp_read);
}





