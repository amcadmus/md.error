#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>
#include "StringSplit.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char * argv[])
{  
  po::options_description desc ("Allow options");

  double xlo, xup;
  unsigned colume;
  std::string file;
  
  desc.add_options()
      ("help,h", "print this message")
      ("colume,c", po::value<unsigned > (&colume)->default_value(2), "colume")
      ("xlow,l", po::value<double > (&xlo)->default_value(5.0), "value of lower bound")
      ("xup,u",  po::value<double > (&xup)->default_value(10.), "value of upper bound")
      ("input,i", po::value<std::string > (&file)->default_value ("force.ref"), "input file");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  int ci = colume - 1;
  double sum = 0.;
  int count = 0;
  char line[20480];
  FILE * fp = fopen (file.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << file << std::endl;
    return 1;
  }
  while (fgets(line, sizeof(line), fp) != NULL){
    if (line[0] == '#') continue;
    std::vector<std::string > words;
    StringOperation::split (std::string(line), words);
    if (words.size() < colume){
      std::cerr << "wrong line fomat, do nothing" << std::endl;
      return 1;
    }
    double xx = atof(words[0].c_str());
    if (xx >= xlo && xx <= xup){
      double value = atof(words[ci].c_str());
      sum += value * value;
      count ++;
    }
  }

  if (count != 0){
    printf ("average_error: %e\n", sqrt(sum / double(count)));
  }

  return 0;
}
