#define FMT_HEADER_ONLY
#include "diagram.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void vertex3::Build(int order) {
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
                 Order - 1, // loopNum
                 5,         // loop index of the first internal K of the vertex
                 0,         // tau index of the InTL leg
                 Chan, RIGHT);
    for (auto &t : Vertex.Tpair) {
      // cout << t[0] << ", " << t[1] << ", " << t[2] << ", " << t[3] << endl;
      int inR = G[0].AddTidxPair({ExtTauIdx, t[INR]});
      int outR = G[1].AddTidxPair({t[OUTR], ExtTauIdx});
      Gidx.push_back(array<int, 2>({inR, outR}));
    }
  }
};

double vertex3::Evaluate() {
  double Factor = 1.0 / pow(2.0 * PI, D);
  // normalization
  if (Order == 0)
    return 1.0;

  // loop order >=1
  vertex4 &Ver4 = Vertex;
  G[0].K = Var.LoopMom[4];
  G[1].K = Var.LoopMom[4];

  for (auto &g : G)
    g.Evaluate();

  Vertex.Evaluate(Var.LoopMom[INL], Var.LoopMom[INR], Var.LoopMom[4],
                  Var.LoopMom[4], false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &gidx = Gidx[i];
    double temp = (SPIN * Ver4.Weight[i][DIR] + Ver4.Weight[i][EX]);
    // attach four G
    for (int j = 0; j < 2; ++j)
      temp *= G[j][gidx[j]];

    Weight += temp;
  }
  return Weight * Factor;
}