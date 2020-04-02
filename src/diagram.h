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
  vector<array<int, 4>> Gidx; // external T list
  // there is only one Tau: from 0 to 1
  ver::weightMatrix Weight;
};

struct sigma {
  int TauNum;
  int LoopNum;
  int OutT;
  ver4 Vertex;
  vector<green> G;
  momentum K; // one internal K is needed
  vector<int> T;
  vector<ver::weightMatrix> Weight;
};

polar BuildPolar(int LoopNum, array<momentum *, 4> ExtLegK);
sigma BuildSigma(int LoopNum, momentum *ExtK);

} // namespace dse
#endif