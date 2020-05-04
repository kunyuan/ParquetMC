#ifndef polar_H
#define polar_H

#include "vertex4.h"

namespace diag {

class polar {
public:
  int Order;
  vertex4 Vertex;
  array<green, 4> G;
  vector<array<int, 4>> Gidx; // external T list

  int TauNum() { return Order + 1; }
  int LoopNum() { return Order; }
  void Build(int Order);
  double Evaluate();
};

class sigma {
  int Order;
  green G1, G2, G3;
  double Weight;
  vector<bubble> Bubble;

  // Ver4, G0, G1, G2 index
  vector<array<int, 4>> Tpair;

  int TauNum() { return Order; }
  int LoopNum() { return Order; }
  void Build(int Order);
  double Evaluate();
};

} // namespace diag
#endif