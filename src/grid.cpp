#include "grid.h"
#include "utility/abort.h"
#include "utility/utility.h"
#include <iomanip>
#include <iostream>
using namespace std;

void logGrid::Initialize(array<double, 2> bound, array<double, 2> idx,
                         double lambda, bool dense2sparse) {
  Bound = bound;
  Idx = idx;
  Lambda = lambda;

  ASSERT_ALLWAYS(bound[0] < bound[1],
                 "Bound[0] should be smaller than Bound[1]!");
  ASSERT_ALLWAYS(idx[0] < idx[1], "Idx[0] should be smaller than Idx[1]!");
  ASSERT_ALLWAYS(Lambda > 0.0, "Lambda must be positive!");

  if (dense2sparse == false) {
    swap(Bound[0], Bound[1]);
    swap(Idx[0], Idx[1]);
    Lambda = -Lambda;
  }

  double l0 = 1.0;
  double l1 = exp(Lambda * (Idx[1] - Idx[0]));
  _b = (Bound[1] - Bound[0]) / (l1 - l0);
  _a = (Bound[0] * l1 - Bound[1] * l0) / (l1 - l0);

  // cout << "Bound: " << Bound[0] << "->" << Bound[1] << endl;
  // cout << "Idx: " << Idx[0] << "->" << Idx[1] << endl;
  // cout << "_a: " << _a << ", _b: " << _b << endl;
  // cout << _a + _b << ", " << _a + _b * exp(-Lambda) << endl;
  // cout << endl;

  ASSERT_ALLWAYS(lambda > 0, "Lambda must be larger than 0!");
};

int logGrid::Floor(double x) {
  ASSERT(x >= Bound[0] && x <= Bound[1], "x must be within the bounds!");
  int idx = Idx[0] + 1.0 / Lambda * log((x - _a) / _b);

  return idx;
}

double logGrid::Grid(int idx) {
  // cout << _a + _b * exp(Lambda * (idx - Idx[0])) << endl;
  return _a + _b * exp(Lambda * (idx - Idx[0]));
}

void tauGrid::Initialize(double Beta, int size, double scale) {

  Size = size;
  Grid.resize(size);
  Weight.resize(size);
  double lambda = Beta / scale / (Size / 2.0);

  _Grid0.Initialize({0.0, Beta / 2.0}, {0.0, Size / 2 - 0.5}, lambda, true);
  _Grid1.Initialize({Beta / 2.0, Beta}, {Size / 2 - 0.5, Size - 1.0}, lambda,
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
  ss << setprecision(8);
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";

  return ss.str();
}

void kFermiGrid::Initialize(double kf, double maxK, int size, double scale) {
  Size = size;
  Grid.resize(size);
  MaxK = maxK;
  kF = kf;
  kFIdx = Size * log(kF) / log(maxK - kF);
  double lambda = kF / scale / kFIdx;

  // the last point of _Grid0 should not be kF!
  _Grid0.Initialize({0.0, kF}, {0.0, kFIdx * 1.0}, lambda, false);
  // the first point of _Grid1 should be kF!
  _Grid1.Initialize({kF, maxK}, {kFIdx * 1.0, Size - 1.0}, lambda, true);

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
  ss << setprecision(8);
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";
  return ss.str();
}

void kBoseGrid::Initialize(double kf, double maxK, int size, double scale) {
  Size = size;
  Grid.resize(size);
  MaxK = maxK;
  kF = kf;
  kFIdx = size / 3;
  TwokFIdx = size / 3 * 2;
  double lambda = kF / scale / kFIdx;

  ASSERT_ALLWAYS(maxK > 2.0 * kF, "MaxK must be larger than the 2kF!");

  _Grid0.Initialize({0.0, kF}, {0.0, kFIdx * 1.0}, lambda, true);
  _Grid1.Initialize({kF, 2.0 * kF}, {kFIdx * 1.0, TwokFIdx * 1.0}, lambda,
                    false);
  _Grid2.Initialize({2.0 * kF, maxK}, {TwokFIdx * 1.0, Size - 1.0}, lambda,
                    true);

  for (int i = 0; i < kFIdx; ++i)
    Grid[i] = _Grid0.Grid(i);

  for (int i = kFIdx; i < TwokFIdx; ++i)
    Grid[i] = _Grid1.Grid(i);

  for (int i = TwokFIdx; i < Size; ++i)
    Grid[i] = _Grid2.Grid(i);

  Grid[0] = 1.0e-6;
  // cout << "Last K: " << Grid[Size - 1] << endl;
  // cout << "K: " << _Grid2.Grid(Size - 1) << ", MaxK: " << maxK << endl;
  // cout << kFIdx << ", " << TwokFIdx << ", " << Grid.size() << endl;
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
  ss << setprecision(8);
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";
  return ss.str();
}

void uniformGrid::Initialize(array<double, 2> bounds, int size) {
  Size = size;
  Grid.resize(size);
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

string uniformGrid::ToString() {
  stringstream ss;
  ss << setprecision(8);
  for (int i = 0; i < Size; ++i)
    ss << Grid[i] << " ";
  return ss.str();
}