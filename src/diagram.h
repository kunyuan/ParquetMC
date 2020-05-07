#ifndef polar_H
#define polar_H

#include "vertex4.h"

namespace diag {

class polar {
public:
  int Order;
  int ExtTauIdx;
  vertex4 Vertex;
  array<green, 4> G;
  vector<array<int, 4>> Gidx; // external T list

  int TauNum() { return Order + 1; }
  int LoopNum() { return Order; }
  void Build(int Order);
  double Evaluate();
};

struct verPair {
  vertex4 LVer;
  vertex4 RVer;
  // map LVerIdx and RVerIdx to VerIdx
  // LVerIdx, G1idx, G2idx, G2idx
  vector<array<int, 4>> Map;
};

class sigma {
public:
  int Order;
  green G1, G2, G3;
  vector<verPair> Bubble;

  int TauNum() { return Order; }
  int LoopNum() { return Order; }
  void Build(int Order);
  double Evaluate();
  string ToString();

  bool Test();
};

} // namespace diag
#endif