#include "grid.h"
#include <algorithm>
#include <assert.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;
using namespace grid;

vector<double> grid::logGrid(const vector<Coeff> &Coeff,
                             const vector<array<int, 2>> &Range) {
  assert(Coeff.size() == Range.size());
  vector<double> grid;
  for (int i = 0; i < (int)Range.size(); ++i)
    for (int j = Range[i][0]; j < Range[i][1]; ++j)
      grid.push_back(Coeff[i].grid(j));
  return grid;
};

void Coeff::init(array<double, 2> _bound, array<double, 2> _idx, double _lambda,
                 bool dense2sparse) {
  bound = _bound;
  assert(bound[0] < bound[1]);
  idx = _idx;
  assert(idx[0] < idx[1]);
  lambda = _lambda;
  assert(lambda > 0.0);

  if (dense2sparse == false) {
    swap(bound[0], bound[1]);
    swap(idx[0], idx[1]);
    lambda = -lambda;
  }

  double l0 = 1.0;
  double l1 = exp(lambda * (idx[1] - idx[0]));
  _b = (bound[1] - bound[0]) / (l1 - l0);
  _a = (bound[0] * l1 - bound[1] * l0) / (l1 - l0);
};

// return the index of a given value
int Coeff::floor(double x) const {
  return idx[0] + 1.0 / lambda * log((x - _a) / _b);
};
// return the Grid with a given index
double Coeff::grid(int i) const {
  return _a + _b * exp(lambda * (i - idx[0]));
};

vector<double> Tau::build(double beta, int _size, double scale) {
  size = _size;
  assert(size > 2);
  assert(size % 2 == 0);
  double lambda = beta / scale / (_size / 2.0);

  _coeff0.init({0.0, beta / 2.0}, {0.0, size / 2 - 0.5}, lambda, true);
  array<int, 2> range0 = {0, size / 2};
  _coeff1.init({beta / 2.0, beta}, {size / 2 - 0.5, size - 1.0}, lambda, false);
  array<int, 2> range1 = {size / 2, size};

  weight.resize(size);
  for (auto &w : weight)
    w = 1.0;

  grid = logGrid({_coeff0, _coeff1}, {range0, range1});

  grid[0] = 1.0e-8;
  grid[size - 1] = beta - 1.0e-8;
  //////   some simple test ////////
  assert(floor(grid[1]) == 1);
  assert(floor((grid[size - 2] + grid[size - 1]) / 2.0) == size - 2);
  assert(floor(0.0) == 0);
  assert(floor(grid[size - 1]) == size - 2);
  //////////////////////////////////
  return grid;
};

int Tau::floor(double x) const {
  // return the index of a given value
  if (x >= grid[1] && x < grid[size / 2 - 1])
    return _coeff0.floor(x);
  else if (x >= grid[size / 2] && x < grid[size - 2])
    return _coeff1.floor(x);
  else if (x >= grid[size - 2])
    return size - 2;
  else if (x>= grid[size/2-1] && x<grid[size/2])
    return size/2-1;
  else
    return 0; // x<Grid[1]
};

string Tau::str() {
  stringstream ss;
  ss << setprecision(12);
  for (auto &g : grid)
    ss << g << " ";
  return ss.str();
};

void FermiK::build(double kF, double maxK, int _size, double scale) {
  assert(maxK > kF);
  size = _size;
  assert(size > 2);
  kFidx = size / 2;
  double lambda = kF / scale / kFidx;

  // the last point of _Grid0 should not be kF!
  _coeff0.init({0.0, kF}, {0.0, kFidx * 1.0}, lambda, false);
  array<int, 2> range0 = {0, kFidx};
  // the first point of _Grid1 should be kF!
  _coeff1.init({kF, maxK}, {kFidx * 1.0, size - 1.0}, lambda, true);
  array<int, 2> range1 = {kFidx, size};

  grid = logGrid({_coeff0, _coeff1}, {range0, range1});

  grid[0] = 1.0e-6;
  // cout << floor(maxK) << endl;
  //////   some simple test ////////
  assert(floor(grid[1]) == 1);
  assert(floor((grid[size - 2] + grid[size - 1]) / 2.0) == size - 2);
  assert(floor(0.0) == 0);
  assert(floor(grid[size - 1]) == size - 2);
  //////////////////////////////////
};
int FermiK::floor(double x) const {
  // return the index of a given value
  if (x >= grid[1] && x < grid[kFidx])
    return _coeff0.floor(x);
  else if (x >= grid[kFidx] && x < grid[size - 2])
    return _coeff1.floor(x);
  else if (x >= grid[size - 2])
    return size - 2;
  else
    return 0; // x<Grid[1]
};
string FermiK::str() {
  stringstream ss;
  ss << setprecision(12);
  for (auto &g : grid)
    ss << g << " ";
  return ss.str();
};

void BoseK::build(double kF, double maxK, int _size, double scale) {
  size = _size;
  kFidx = size / 3;
  twokFidx = size / 3 * 2;
  double lambda = kF / scale / kFidx;
  assert(size > 3);
  assert(maxK > 2.0 * kF);

  _coeff0.init({0.0, kF}, {0.0, kFidx * 1.0}, lambda, true);
  array<int, 2> range0 = {0, kFidx};
  _coeff1.init({kF, 2.0 * kF}, {kFidx * 1.0, twokFidx * 1.0}, lambda, false);
  array<int, 2> range1 = {kFidx, twokFidx};
  _coeff2.init({2.0 * kF, maxK}, {twokFidx * 1.0, size - 1.0}, lambda, true);
  array<int, 2> range2 = {twokFidx, size};

  grid = logGrid({_coeff0, _coeff1, _coeff2}, {range0, range1, range2});

  grid[0] = 1.0e-6;
  //////   some simple test ////////
  assert(floor(grid[1]) == 1);
  assert(floor((grid[size - 2] + grid[size - 1]) / 2.0) == size - 2);
  assert(floor(0.0) == 0);
  assert(floor(grid[size - 1]) == size - 2);
  //////////////////////////////////
};
int BoseK::floor(double x) const {
  // return the index of a given value
  if (x >= grid[1] && x < grid[kFidx])
    return _coeff0.floor(x);
  else if (x >= grid[kFidx] && x < grid[twokFidx])
    return _coeff1.floor(x);
  else if (x >= grid[twokFidx] && x < grid[size - 2])
    return _coeff2.floor(x);
  else if (x >= grid[size - 2])
    return size - 2;
  else
    return 0; // x<Grid[1]
};
string BoseK::str() {
  stringstream ss;
  ss << setprecision(12);
  for (auto &g : grid)
    ss << g << " ";
  return ss.str();
};

void Uniform::build(std::array<double, 2> bounds, int _size) {
  size = _size;
  assert(size > 1);
  grid.resize(size);
  lowerBound = bounds[0];
  delta = (bounds[1] - bounds[0]) / (size - 1);
  for (int i = 0; i < size; ++i)
    grid[i] = i * delta + lowerBound;

  assert(floor(grid[size - 1]) == size - 2);
  assert(floor(grid[0]) == 0);
};
int Uniform::floor(double x) const {
  if (x <= grid[size - 2])
    return (x - lowerBound ) / delta;
  else
    return size - 2; // x>grid[size-2]
};
std::string Uniform::str() {
  stringstream ss;
  ss << setprecision(12);
  for (auto &g : grid)
    ss << g << " ";
  return ss.str();
};
