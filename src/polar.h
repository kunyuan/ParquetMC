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

} // namespace diag
#endif