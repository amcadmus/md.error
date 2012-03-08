#ifndef __VectorOperation_h_wanghan__
#define __VectorOperation_h_wanghan__

#include <vector>
#include <string>
#include <cmath>

namespace VectorOperation{

// level 1
    // set all v component 0
    inline void zero (std::vector<double > & v)
    {
      std::fill (v.begin(), v.end(), 0.);
    }
    // v = alpha * v
    inline void scale (std::vector<double > & v, const double & alpha)
    {
      v[0]*=alpha;
      v[1]*=alpha;
      v[2]*=alpha;
    }
    // v += alpha * a
    inline void add (std::vector<double > & v, 
	      const double & alpha, const std::vector<double > & a)
    {
      v[0] += alpha * a[0];
      v[1] += alpha * a[1];
      v[2] += alpha * a[2];
    }
    // v += alpha * a + beta * b
    inline void add (std::vector<double > & v,
	      const double & alpha, const std::vector<double > & a,
	      const double & beta,  const std::vector<double > & b)
    {
      v[0] += alpha * a[0] + beta * b[0];
      v[1] += alpha * a[1] + beta * b[1];
      v[2] += alpha * a[2] + beta * b[2];
    }
    // v = s * v + alpha * a
    inline void sadd (std::vector<double > & v, const double & s, 
	       const double & alpha, const std::vector<double > & a)
    {
      v[0] = s * v[0] + alpha * a[0];
      v[1] = s * v[1] + alpha * a[1];
      v[2] = s * v[2] + alpha * a[2];
    }
    // v = s * v + alpha * a + beta * b
    inline void sadd (std::vector<double > & v, const double & s,
	       const double & alpha, const std::vector<double > & a,
	       const double & beta,  const std::vector<double > & b)
    {
      v[0] = s * v[0] + alpha * a[0] + beta * b[0];
      v[1] = s * v[1] + alpha * a[1] + beta * b[1];
      v[2] = s * v[2] + alpha * a[2] + beta * b[2];
    }
    // return u .* v
    inline double dot (const std::vector<double > & u, const std::vector<double > & v)
    {
      return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }
    // print to stdout
    void print (const std::vector<double > & v, 
		const bool & inColume = true);
    // print to file
    void print (const std::string & filename, const std::vector<double > & v);
    // read from a file
    bool read (const std::string & filename, std::vector<double > & v);


};



#endif

