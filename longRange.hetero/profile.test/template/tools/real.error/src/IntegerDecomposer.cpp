#include "ToolBox.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>

struct vector_head_less_than 
{
  bool operator () (const std::vector<int > & a,
		    const std::vector<int > & b) {
    return a[0] < b[0];
  }
};

int intpow (int a, int index)
{
  if (index <= 0){
    return 1;
  }
  else {
    int re = 1;
    while (index -- >= 1){
      re *= a;
    }
    return re;
  }
}

ToolBox::IntegerDecomposer::IntegerDecomposer ()
{
  good.clear();
  for (int a = 0; a <= int(log(256)/log(2)); ++ a){
    for (int b = 0; b <= int(log(256)/log(3)); ++ b){
      for (int c = 0; c <= int(log(256)/log(5)); ++ c){
	for (int d = 0; d <= int(log(256)/log(7)); ++ d){
	  int tmp = intpow(2, a) * intpow(3, b) * intpow(5, c) * intpow(7, d);
	  if (tmp >= 1 && tmp <= 255){
	    std::vector<int > data;
	    data.push_back(tmp);
	    data.push_back(a);
	    data.push_back(b);
	    data.push_back(c);
	    data.push_back(d);
	    good.push_back(data);
	  }
	}
      }
    }
  }
  vector_head_less_than less;
  std::sort (good.begin(), good.end(), less);
//   std::vector<std::vector<int > >::iterator i = good.begin();
//   for (; i != good.end(); ++i){
//     std::copy (i->begin(), i->end(), std::ostream_iterator<int >(std::cout, "\t"));
//     std::cout << std::endl;
//   }
}

int ToolBox::IntegerDecomposer::decompose ( 
    int input, std::vector<int > & result)
{
  result.resize(5, 0);
  result[0] = 1;
  
 // double n16;
  if (input >= 256) {
    double n16d = log(input) / log(16);
    int n16i = ((n16d == int(n16d)) ? int(n16d) : int(n16d) + 1);
    result[1] += 4 * (n16i - 2);
    result[0] *= intpow(16, (n16i - 2)); 
    input = int (double(input) / double(result[0]) + 0.5);
  }

  std::vector<std::vector<int > >::iterator select;
  std::vector<std::vector<int > >::iterator i0 = good.begin();
  std::vector<std::vector<int > >::iterator i1 = good.begin();
  i1 ++;
  for (; i1 != good.end(); ++i1, ++i0){
    if (i1->front() >= input && i0->front() <= input){
      fabs (i1->front() - input) <= fabs (i0->front() - input) ?
	  select = i1 :
	  select = i0;
      result[0] *= select->operator[](0);
      result[1] += select->operator[](1);
      result[2] += select->operator[](2);
      result[3] += select->operator[](3);
      result[4] += select->operator[](4);
      return result[0];
    }
  }
  return -1;
}

int ToolBox::IntegerDecomposer::decompose (int input)
{
  std::vector<int > result (5, 0);
  result[0] = 1;
  
//   double n16;
  if (input >= 256) {
    double n16d = log(input) / log(16);
    int n16i = ((n16d == int(n16d)) ? int(n16d) : int(n16d) + 1);
    result[1] += 4 * (n16i - 2);
    result[0] *= intpow(16, (n16i - 2)); 
    input = int (double(input) / double(result[0]) + 0.5);
  }

  std::vector<std::vector<int > >::iterator select;
  std::vector<std::vector<int > >::iterator i0 = good.begin();
  std::vector<std::vector<int > >::iterator i1 = good.begin();
  i1 ++;
  for (; i1 != good.end(); ++i1, ++i0){
    if (i1->front() >= input && i0->front() <= input){
      fabs (i1->front() - input) <= fabs (i0->front() - input) ?
	  select = i1 :
	  select = i0;
      result[0] *= select->operator[](0);
      result[1] += select->operator[](1);
      result[2] += select->operator[](2);
      result[3] += select->operator[](3);
      result[4] += select->operator[](4);
      return result[0];
    }
  }
  return -1;
}

   
