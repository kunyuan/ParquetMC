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
  ver4 Vertex;
  momentum KOutL, KOutR;
  array<vector<green>, 4> G;
  vector<array<int, 4>> Gidx; // external T list
  // there is only one Tau: from 0 to 1
  double Weight;
};

struct sigma {
  int TauNum;
  int LoopNum;
  ver4 Vertex;
  momentum *K; // one internal K is needed
  vector<green> G;

  vector<int> T;
  vector<double> Weight;

  // map Vertex weight idx to G idx and sigma idx
  vector<int> Gidx;
  vector<int> SigTidx;
};

polar BuildPolar(int LoopNum, momentum *InL, momentum *InR);
sigma BuildSigma(int LoopNum, momentum *ExtK, momentum *InterK);

void SetSigmaMom(sigma &Sigma, momentum *ExtK, momentum *InterK);

} // namespace dse
#endif