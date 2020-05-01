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
  // vector<green> GInL, GOutL, GInR, GOutR;
  array<vector<green>, 4> G;
  vector<array<int, 4>> Gidx; // external T list
  // there is only one Tau: from 0 to 1
  double Weight;
};

struct sigma {
  int TauNum;
  int LoopNum;
  ver4 Vertex;
  vector<green> G;

  vector<int> T;
  vector<double> Weight;

  // map Vertex weight idx to G idx and sigma idx
  vector<int> Gidx;
  vector<int> SigTidx;
};

polar BuildPolar(int LoopNum);
sigma BuildSigma(int LoopNum);

} // namespace dse
#endif