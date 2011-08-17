#include "VectorOperation.h"
#include <iostream>
#include <cassert>
#include <fstream>

void VectorOperation::print (const std::vector<double > & v, const bool & inColume)
{
  if (inColume){
    for (std::vector<double >::const_iterator pv = v.begin(); pv != v.end();){
      std::cout << *(pv ++) << std::endl;
    }
  }
  else{
    for (std::vector<double >::const_iterator pv = v.begin(); pv != v.end();){
      std::cout << *(pv ++) << "\t";
    }
    std::cout << std::endl;
  }
}

void VectorOperation::print (const std::string & filename, const std::vector<double > & v)
{
  FILE * fp = fopen (filename.c_str(), "w");
  std::vector<double >::const_iterator it = v.begin();
  for (; it != v.end();){
    fprintf(fp, "%.16e\n", *(it++));
  }
  fclose (fp);
} 

bool VectorOperation::read (const std::string & filename, std::vector<double > & v)
{
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL) return false;
  double tmp;
  v.clear ();
  int info;
  
  while ((info = fscanf (fp, "%lf", &tmp)) != EOF){
    assert (info == 1);
    v.push_back (tmp);
  }
  
  fclose (fp);
  return true;
}
