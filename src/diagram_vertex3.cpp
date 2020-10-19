#define FMT_HEADER_ONLY
#include "diagram.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void polar::Build(int order) {
  ASSERT_ALLWAYS(order >= 0, "Polar order must be larger than 0!");
  Order = order;
  ExtTauIdx = MaxTauNum - 1;
  ASSERT_ALLWAYS(ExtTauIdx > TauNum(), "MaxTauNum is too small!");

  // vertex is only needed for order>=2
  if (Order >= 2) {
    vector<channel> Chan = {I, T, U, S, TC, UC};
    // vector<channel> Chan = {
    //     T,
    // };
    Vertex.Build(0,         // level
                 Order - 2, // loopNum
                 3,         // loop index of the first internal K of the vertex
                 1,         // tau index of the InTL leg
                 Chan, RIGHT);
    for (auto &t : Vertex.Tpair) {
      // cout << t[0] << ", " << t[1] << ", " << t[2] << ", " << t[3] << endl;
      int inR = G[0].AddTidxPair({ExtTauIdx, t[INR]});
      int outR = G[1].AddTidxPair({t[OUTR], ExtTauIdx});
      Gidx.push_back(array<int, 4>({inR, outR}));
    }
  }
};

double polar::Evaluate() {
  double Factor = 1.0 / pow(2.0 * PI, D);
  // normalization
  if (Order == 0)
    return 1.0;
  else if (Order == 1) {
    return 1.0;
  }

  // loop order >=2
  vertex4 &Ver4 = Vertex;
  G[INL].K = Var.LoopMom[1];
  G[INR].K = Var.LoopMom[2];
  G[OUTL].K = Var.LoopMom[1] - Var.LoopMom[0];
  G[OUTR].K = Var.LoopMom[2] + Var.LoopMom[0];

  for (auto &g : G)
    g.Evaluate();

  Vertex.Evaluate(G[INL].K, G[OUTL].K, G[INR].K, G[OUTR].K, false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &gidx = Gidx[i];
    double temp =
        (SPIN * SPIN * Ver4.Weight[i][DIR] + SPIN * Ver4.Weight[i][EX]);
    // attach four G
    for (int j = 0; j < 4; ++j)
      temp *= G[j][gidx[j]];

    Weight += temp;
  }
  return Weight * Factor * Factor;
}