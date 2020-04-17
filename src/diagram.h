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
  momentum KOutL, KOutR;
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
  momentum *K; // one internal K is needed
  vector<green> G;

  vector<int> T;
  vector<double> Weight;

  // map Vertex Tidx and G idx to sigma Tidx
  vector<int> Gidx;
  vector<int> VerTidx;
  vector<int> SigTidx;
};

polar BuildPolar(int LoopNum, momentum *InL, momentum *InR);
sigma BuildSigma(int LoopNum, momentum *ExtK, momentum *InterK);

} // namespace dse
#endif