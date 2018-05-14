#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

int main(int argc, char **argv) {
  std::ofstream ofs("./test_cpp.ppm", std::ios::out | std::ios::binary);
  // ofs << (unsigned char)(float(2.0));
  ofs << "P6\n" << 640 << " " << 480 << "\n255\n";
}
