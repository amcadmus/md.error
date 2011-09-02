#include "AssignFcorr.h"
#include <fstream>


AssignFcorr::
AssignFcorr ()
    : malloced (false)
{
}

AssignFcorr::
~AssignFcorr ()
{
  freeAll();
}


void AssignFcorr::
freeAll ()
{
  if (malloced) {
    free (hfx);
    free (hfy);
    free (hfz);
    malloced = false;
  }
}



void AssignFcorr::    
print_xy (const char * file) const 
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    for (int k = 0; k < nz; ++k){
      fprintf (fp, "%f %f  %e   %e %e %e\n",
	       (i + 0.5) * box[0] / double(nx),
	       (k + 0.5) * box[2] / double(nz),
	       sqrt(
		   hfx [index3to1(i, ny/2, k)] * hfx [index3to1(i, ny/2, k)] +
		   hfy [index3to1(i, ny/2, k)] * hfy [index3to1(i, ny/2, k)] +
		   hfz [index3to1(i, ny/2, k)] * hfz [index3to1(i, ny/2, k)] ),
	       hfx [index3to1(i, ny/2, k)],
	       hfy [index3to1(i, ny/2, k)],
	       hfz [index3to1(i, ny/2, k)]
	  );
    }
    fprintf (fp, "\n");
  }
  fclose (fp);
}

void AssignFcorr::
init_read (const char * file)
{
  fp_read = fopen (file, "r");
  if (fp_read == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit(1);
  }
  int tmpnn[3];
  int status;
  
  status = fread (box, sizeof(double), 3, fp_read);
  if (status != 3){
    printf ("format error\n");
    exit(1);
  }
  status = fread (tmpnn, sizeof(int), 3, fp_read);
  if (status != 3){
    printf ("format error\n");
    exit(1);
  }

  nx = tmpnn[0];
  ny = tmpnn[1];
  nz = tmpnn[2];
  nele = nx * ny * nz;

  freeAll();
  hfx = (float *) malloc (sizeof(float ) * nele);
  hfy = (float *) malloc (sizeof(float ) * nele);
  hfz = (float *) malloc (sizeof(float ) * nele);
  malloced = true;
}

bool AssignFcorr::
read (float & time) 
{
  int status;
  status = fread (&time, sizeof(float), 1, fp_read);
  if (status != 1){
    printf ("format error or reached end\n");
    return false;
  }
  status = fread (hfx, sizeof(float), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    return false;
  }
  status = fread (hfy, sizeof(float), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    return false;
  }
  status = fread (hfz, sizeof(float), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    return false;
  }
  return true;
}

void AssignFcorr::
end_read () const
{
  fclose (fp_read);
}

