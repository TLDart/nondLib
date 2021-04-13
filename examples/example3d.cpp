#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

#include "../nondlib.hpp"

template <typename T, typename RngGen>
std::vector<std::vector<T>> generatePointset(int sz, int dim, RngGen& gen) {
  std::uniform_int_distribution<T> dis(0, 10000);

  std::vector<std::vector<T>> pointset;
  for (int i = 0; i < sz; i++) {
    std::vector<T> row;
    for (int j = 0; j < dim; j++) {
      row.push_back(dis(gen));
    }
    pointset.push_back(std::move(row));
  }
  return pointset;
}

template <typename T>
void printPoints(T const& points) {
  for (auto const& p : points) {
    for (auto c : p) {
      std::cout << c << " ";
    }
    std::cout << "\n";
  }
}

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());

  auto set = generatePointset<int>(100, 3, gen);
  auto set2 = set;
  auto set3 = set;
  auto set4 = set;
  auto set5 = set;
  auto mx = std::vector<int>{1, 1, 1};

  std::cout << "# Before filter\n";
  printPoints(set);

  nondlib::inplace::filterDimSweep3D(set, mx, 0);
  std::sort(set.begin(), set.end());
  std::cout << "\n# After filter\n";
  printPoints(set);

  auto set6 = nondlib::notinplace::filterDimSweep3D(set2, mx, 0);
  auto set7 = nondlib::notinplace::filterDimSweep3D(set2, mx, 1);
  auto set8 = nondlib::notinplace::filterDimSweep3D(set2, mx, 2);
  auto set9 = nondlib::notinplace::filterQuadD(set2, mx);
  auto set10 = nondlib::notinplace::filterDivConqDC(set2, mx);
  nondlib::inplace::filterQuadD(set2, mx);
  nondlib::inplace::filterDivConqDCC(set3, mx);
  nondlib::inplace::filterDimSweep3D(set4, mx, 1);
  nondlib::inplace::filterDimSweep3D(set5, mx, 2);

  std::sort(set2.begin(), set2.end());
  std::sort(set3.begin(), set3.end());
  std::sort(set4.begin(), set4.end());
  std::sort(set5.begin(), set5.end());
  std::sort(set6.begin(), set6.end());
  std::sort(set6.begin(), set6.end());
  std::sort(set7.begin(), set7.end());
  std::sort(set8.begin(), set8.end());
  std::sort(set9.begin(), set9.end());
  std::sort(set10.begin(), set10.end());

  assert(set == set2);
  assert(set == set3);
  assert(set == set4);
  assert(set == set5);
  assert(set == set6);
  assert(set == set7);
  assert(set == set8);
  assert(set == set9);
  assert(set == set10);

  auto newpoints = generatePointset<int>(100, 3, gen);
  for (auto& p : newpoints) {
    set = nondlib::notinplace::updateMaximaND(set, mx, p);
    nondlib::inplace::updateMaximaND(set2, mx, p);
  }
  assert(set == set2);

  std::cout << "\n# After update\n";
  printPoints(set);

  return 0;
}
