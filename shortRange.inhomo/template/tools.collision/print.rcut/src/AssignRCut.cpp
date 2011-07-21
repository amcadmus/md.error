#include "AssignRCut.h"
#include <fstream>


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
uniform (const double & rc)
{
  for (int i = 0; i < nele; ++i){
    hrcut[i] = rc;
  }
}

void AssignRCut::
freeAll ()
{
  if (malloced) {
    free (hrcut);
    malloced = false;
  }
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
	     (i + 0.5) * box[0] / double(nx),
	     hrcut [index3to1(i, 0, 0)]
	);
  }
  fclose (fp);
}


void AssignRCut::    
print_xy (const char * file) const 
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    exit(1);
  }

  for (int i = 0; i < nx; ++i){
    for (int k = 0; k < nz; ++k){
      fprintf (fp, "%f %f %e\n",
	       (i + 0.5) * box[0] / double(nx),
	       (k + 0.5) * box[2] / double(nz),
	       hrcut [index3to1(i, ny/2, k)]
	  );
    }
    fprintf (fp, "\n");
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
  tmpbox[0] = box[0];
  tmpbox[1] = box[1];
  tmpbox[2] = box[2];
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
write (const float & time) const
{
  float tmptime = time;
  fwrite (&tmptime, sizeof(float), 1,    fp_write);
  fwrite (hrcut,    sizeof(float), nele, fp_write);
}

void AssignRCut::
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
  hrcut = (float *) malloc (sizeof(float ) * nele);
  malloced = true;
}

bool AssignRCut::
read (float & time) 
{
  int status;
  status = fread (&time, sizeof(float), 1, fp_read);
  if (status != 1){
    printf ("format error or reached end\n");
    return false;
  }
  status = fread (hrcut, sizeof(float), nele, fp_read);
  if (status != nele){
    printf ("format error\n");
    return false;
  }
  return true;
}

void AssignRCut::
end_read () const
{
  fclose (fp_read);
}

