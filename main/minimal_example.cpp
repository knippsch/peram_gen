#include <complex> 
#include <iostream>
#include <random>

#include "ranlxs.h"

int main(int argc, char *argv[]){

  const size_t length = 5;
  const int seed = 1227;

  rlxs_init(0, seed);
  std::vector<float> rnd(length);
  ranlxs(&(rnd[0]), length);
  std::cout << "Lüscher rnd:\n" << std::endl;
  for(const auto& s : rnd)
    std::cout << s << std::endl;

  std::cout << "\n\nLüscher rnd from std:\n" << std::endl;
  std::ranlux24 g1(seed);
  std::uniform_real_distribution<float> dr(0.0, 1.0);
  std::vector<float> rnd2(length);
  for(auto& s : rnd2)
    std::cout << (s = dr(g1)) << std::endl;

}
