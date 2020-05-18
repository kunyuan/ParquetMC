#include "grid.h"
#include "utility/abort.h"
#include "utility/utility.h"
#include <iostream>
using namespace std;

void logGrid::Initialize(array<double, 2> bound, array<double, 2> idx,
                         double lambda, bool dense2sparse) {
  Bound = bound;
  Idx = idx;
  Lambda = lambda / (abs(idx[1] - idx[0]));

  ASSERT_ALLWAYS(bound[0] < bound[1],
                 "Bound[0] should be smaller than Bound[1]!");
  ASSERT_ALLWAYS(idx[0] < idx[1], "Idx[0] should be smaller than Idx[1]!");
  ASSERT_ALLWAYS(Lambda > 0.0, "Lambda must be positive!");

  if (dense2sparse == false) {
    swap(Bound[0], Bound[1]);
    swap(Idx[0], Idx[1]);
    Lambda = -Lambda;
  }

  double l0 = exp(Lambda * Idx[0]);
  double l1 = exp(Lambda * Idx[1]);
  _b = (Bound[1] - Bound[0]) / (l1 - l0);
  _a = (Bound[0] * l1 - Bound[1] * l0) / (l1 - l0);

  ASSERT_ALLWAYS(lambda > 0, "Lambda must be larger than 0!");
};

int logGrid::Floor(double x) {
  ASSERT(x >= Bound[0] && x <= Bound[1], "x must be within the bounds!");
  int idx = Idx[0] + 1.0 / Lambda * log((x - _a) / _b);

  return idx;
}

double logGrid::Grid(int idx) { return _a + _b * exp(Lambda * (idx - Idx[0])); }

void tauGrid::Initialize(double Beta, int size, double lambda) {

  Size = size;

  _Grid0.Initialize({0.0, Beta / 2.0}, {0.0, Size / 2 - 0.5}, lambda, true);
  _Grid1.Initialize({Beta / 2.0, Beta}, {Size / 2 - 0.5, Size - 1}, lambda,
                    false);

  for (int i = 0; i < Size / 2; ++i) {
    Grid[i] = _Grid0.Grid(i);
    Weight[i] = 1.0;
  }

  for (int i = Size / 2; i < Size; ++i) {
    Grid[i] = _Grid1.Grid(i);
    Weight[i] = 1.0;
  }

  Grid[0] = 1.0e-8;
  Grid[Size - 1] = Beta - 1.0e-8;

  return;
}

int tauGrid::Floor(double x) {
  int idx;

  if (x >= Grid[1] && x < Grid[Size / 2 - 1])
    idx = _Grid0.Floor(x);
  else if (x >= Grid[Size / 2] && x < Grid[Size - 2])
    idx = _Grid1.Floor(x);
  else if (x < Grid[1])
    idx = 0;
  else
    idx = Size - 2; // x>=Grid[Size-2]

  ASSERT(idx < Size && idx >= 0, "Out of range!");
  return idx;
}

string tauGrid::ToString() {
  stringstream ss;
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";

  return ss.str();
}

void kFermiGrid::Initialize(double kf, double maxK, int size, double lambda) {
  Size = size;
  MaxK = maxK;
  kF = kf;
  kFIdx = Size / 3;

  // the last point of _Grid0 should not be kF!
  _Grid0.Initialize({0.0, kF}, {0.0, kFIdx}, lambda, false);
  // the first point of _Grid1 should be kF!
  _Grid1.Initialize({kF, maxK}, {kFIdx, Size - 1}, lambda, true);

  for (int i = 0; i < kFIdx; ++i)
    Grid[i] = _Grid0.Grid(i);

  for (int i = kFIdx; i < Size; ++i)
    Grid[i] = _Grid1.Grid(i);

  Grid[0] = 1.0e-6;
}

int kFermiGrid::Floor(double x) {
  int idx;
  if (x >= Grid[1] && x < Grid[kFIdx])
    idx = _Grid0.Floor(x);
  else if (x >= Grid[kFIdx])
    idx = _Grid1.Floor(x);
  else
    idx = 0; // x<Grid[1]
  ASSERT(idx < Size && idx >= 0, "Out of range!");
  return idx;
}

string kFermiGrid::ToString() {
  stringstream ss;
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";
  return ss.str();
}

void kBoseGrid::Initialize(double kf, double maxK, int size, double lambda) {
  Size = size;
  MaxK = maxK;
  kF = kf;
  lambda = lambda;
  kFIdx = size / 3;
  TwokFIdx = size / 3 * 2;

  _Grid0.Initialize({0.0, kF}, {0.0, kFIdx}, lambda, false);
  _Grid1.Initialize({kF, 2.0 * kF}, {kFIdx, TwokFIdx}, lambda, true);
  _Grid2.Initialize({2.0 * kF, maxK}, {TwokFIdx, Size - 1}, lambda, false);

  for (int i = 0; i < kFIdx; ++i)
    Grid[i] = _Grid0.Grid(i);

  for (int i = kFIdx; i < TwokFIdx; ++i)
    Grid[i] = _Grid1.Grid(i);

  for (int i = TwokFIdx; i < Size - 1; ++i)
    Grid[i] = _Grid1.Grid(i);

  Grid[0] = 1.0e-6;
}

int kBoseGrid::Floor(double x) {
  int idx;
  if (x >= Grid[1] && x < Grid[kFIdx])
    idx = _Grid0.Floor(x);
  else if (x >= Grid[kFIdx] && x < Grid[TwokFIdx])
    idx = _Grid1.Floor(x);
  else if (x >= Grid[TwokFIdx])
    idx = _Grid2.Floor(x);
  else
    idx = 0; // x<Grid[1]
  ASSERT(idx < Size && idx >= 0, "Out of range!");
  return idx;
}

string kBoseGrid::ToString() {
  stringstream ss;
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";
  return ss.str();
}

void uniformGrid::Initialize(array<double, 2> bounds, int size) {
  Size = size;
  LowerBound = bounds[0];
  Delta = (bounds[1] - bounds[0]) / (size - 1);
  for (int i = 0; i < size; ++i)
    Grid[i] = i * Delta;
}

int uniformGrid::Floor(double x) {
  int idx = (LowerBound - x) / Delta;
  ASSERT(idx < Size, "Out of range!");
  return idx;
}