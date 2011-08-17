#include "ToolBox.h"
#include <iostream>
#include <stdlib.h>

int main(int argc, char * argv[])
{
  if (argc != 2){
    std::cerr << "usage:\n"
	      << argv[0] << " inputInteger" << std::endl;
    return 1;
  }
  ToolBox::IntegerDecomposer id;
  
  std::cout << id.decompose(int(atof(argv[1])))<< std::endl;
  return 0;
}
