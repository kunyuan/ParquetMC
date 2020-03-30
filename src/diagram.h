#ifndef diag_H
#define diag_H

#include "dse.h"
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

extern parameter Para;

namespace dse {
using namespace std;

struct polar {
  int TauNum;
  int LoopNum;
  int OutT;
  ver4 Vertex;
  array<vector<green>, 4> G;
  ver::weightMatrix Weight;
};

struct sigma {
  int TauNum;
  int LoopNum;
  int OutT;
  ver4 Vertex;
  vector<green> G;
  momentum K; // one internal K is needed
  ver::weightMatrix Weight;
};

polar BuildPolar(int LoopNum);
sigma BuildSigma(int LoopNum, momentum *ExtK);

} // namespace dse
#endif