#ifndef __wanghan_ToolBox_h__
#define __wanghan_ToolBox_h__
#include<vector>

typedef double value_type;

namespace ToolBox{    
    void init_genrand(unsigned long s);
    unsigned long genrand_int32(void);
    long genrand_int31(void);
    double genrand_real1(void); // in [0,1]
    double genrand_real2(void); // in [0,1)
    double genrand_real3(void); // in (0,1)
    double genrand_res53(void);
    
    struct Integral3DInfo;
    
    template <typename TrinaryFunction, typename ValueType >
    class Integral3D;

    class IntegerDecomposer;
}

class ToolBox::IntegerDecomposer
{
  std::vector<std::vector<int > > good; 
public:
  IntegerDecomposer ();
  int decompose (int input, std::vector<int > & index);
  int decompose (int input);
}
    ;



#endif
