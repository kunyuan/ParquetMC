#include "grid.h"
#include "global.h"
#include <iostream>
using namespace std;

void grid::Initialize(array<double, 2> bound, int size, bool dense2Sparse,
                      double lambda) {
  Bound = bound;
  Size = size;
  Lambda = lambda;
  Dense2Sparse = dense2Sparse;

  //   Gamma = log(1.0 + exp(lambda) * (Bound[1] / Bound[0] - 1.0));
  _Gamma = (Bound[1] - Bound[0]) / (exp(Lambda) - 1.0);
  _Factor = Size / Lambda;

  ASSERT_ALLWAYS(bound[0] < bound[1],
                 "Bound[0] should be smaller than Bound[1]!");
  ASSERT_ALLWAYS(lambda > 0, "Lambda must be larger than 0!");

  double grid;
  for (int i = 0; i < Size; ++i) {

    if (Dense2Sparse)
      grid = Bound[0] + _Gamma * (exp(Lambda * i / Size) - 1.0);
    else
      grid = Bound[1] + _Gamma * (1.0 - exp(Lambda * (Size - i) / Size));

    _Grid.push_back(grid);
  }
};

array<int, 2> grid::Index(double x) {
  ASSERT(x >= Bound[0] && x <= Bound[1], "x must be within the bounds!");
  int idx;
  if (Dense2Sparse)
    idx = _Factor * log(1.0 + (x - Bound[0]) / _Gamma);
  else
    idx = _Factor * log(1.0 + (Bound[1] - x) / _Gamma);
  return array<int, 2>({idx, idx + 1});
}

double &grid::Grid(int index) { return _Grid[index]; };

std::string grid::ToString() {
  stringstream ss;
  for (auto &g : _Grid)
    ss << g << " ";
  return ss.str();
};

void tauGrid::Initialize(double Beta, int size, double lambda) {
  Size = size;
  _Grid[0].Initialize({0.0, Beta / 2.0}, Size / 2 + 1, true, lambda);
  _Grid[1].Initialize({Beta / 2.0, Beta}, Size / 2, false, lambda);
  _Grid[0].Grid(0) = 1.0e-8;
  _Grid[1].Grid(Size / 2 - 1) = Beta - 1.0e-8;
  return;
}

double tauGrid::Grid(int Index) {
  ASSERT(Index < Size, "Out of range!");
  if (Index < Size / 2)
    return _Grid[0].Grid(Index);
  else
    return _Grid[1].Grid(Index - Size / 2);
}

string tauGrid::ToString() {
  stringstream ss;
  for (int i = 0; i < Size; ++i)
    ss << Grid(i) << " ";

  return ss.str();
}
