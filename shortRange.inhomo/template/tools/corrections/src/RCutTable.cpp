#include "RCutTable.h"

RCutTable::
RCutTable ()
{
}

RCutTable::
~RCutTable ()
{
}


void RCutTable::
load_rc (const char * file)
{
  FILE * fp = fopen (file, "r");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit (1);
  }

  fscanf (fp, "%d", &nrc);
  rcList.resize (nrc);
  for (int i = 0; i < nrc; ++i){
    fscanf (fp, "%lf", &rcList[i]);
  }
  boxsize.resize(3);
  fscanf (fp, "%lf%lf%lf", &boxsize[0], &boxsize[1], &boxsize[2]);
  fscanf (fp, "%d%d%d", &nx, &ny, &nz);
  hx = boxsize[0] / nx;
  hy = boxsize[1] / ny;
  hz = boxsize[2] / nz;
  nele = nx * ny * nz;
  rcutIndex.resize (nele);
  for (int i = 0; i < nele; ++i){
    fscanf (fp, "%d", &rcutIndex[i]);
  }
  
  // freeAll ();
  // fread (&nrc, sizeof(int), 1, fp);
  // rcList = (double *) malloc (sizeof(double) * nrc);
  // fread (rcList, sizeof(double), nrc, fp);
  // fread (boxsize, sizeof(double), 3, fp);
  // fread (&nx, sizeof(int), 1, fp);
  // fread (&ny, sizeof(int), 1, fp);
  // fread (&nz, sizeof(int), 1, fp);
  // hx = boxsize[0] / nx;
  // hy = boxsize[1] / ny;
  // hz = boxsize[2] / nz;
  // nele = nx * ny * nz;
  // rcutIndex = (int *) malloc (sizeof(int) * nele);
  // fread (rcutIndex, sizeof(int), nele, fp);
  // malloced = true;
  
  fclose (fp);
}

void RCutTable::
save_rc (const char * file) const
{
  FILE * fp = fopen (file, "w");
  if (fp == NULL){
    fprintf (stderr, "cannot open file %s\n", file);
    exit (1);
  }

  fprintf (fp, "%d ", nrc);
  for (int i = 0; i < nrc; ++i){
    fprintf (fp, "%f ", rcList[i]);
  }
  fprintf (fp, "\n%f %f %f\n", boxsize[0], boxsize[1], boxsize[2]);
  fprintf (fp, "%d %d %d\n", nx, ny, nz);
  for (int i = 0; i < nele; ++i){
    fprintf (fp, "%d ", rcutIndex[i]);
  }
  fprintf(fp, "\n");
  
  // fwrite (&nrc, sizeof(int), 1, fp);
  // fwrite (rcList, sizeof(double), nrc, fp);
  // fwrite (boxsize, sizeof(double), 3, fp);
  // fwrite (&nx, sizeof(int), 1, fp);
  // fwrite (&ny, sizeof(int), 1, fp);
  // fwrite (&nz, sizeof(int), 1, fp);
  // fwrite (rcutIndex, sizeof(int), nele, fp);
  
  fclose (fp);
}

