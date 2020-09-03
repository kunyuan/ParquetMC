#define FMT_HEADER_ONLY
#include "diagram.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void sigma2::Build(int order) {
  ASSERT_ALLWAYS(order >= 0, "Sigma LoopNum must be larger than 0!");
  Order = order;
  ExtTauIdx = MaxTauNum - 1;
  if (order <= 1)
    return;
  ASSERT_ALLWAYS(ExtTauIdx > TauNum(), "MaxTauNum is too small!");

  // vertex is only needed for order>=2
  if (Order >= 2) {
    vector<channel> Chan = {I, T, U, S, TC, UC};
    // vector<channel> Chan = {T, U, S};
    // vector<channel> Chan = {
    //     T,
    // };
    Vertex.Build(0,         // level
                 Order - 2, // loopNum
                 3,         // loop index of the first internal K of the vertex
                 0,         // tau index of the InTL leg
                 Chan, RIGHT);

    for (auto &t : Vertex.Tpair) {
      int t1 = G1.AddTidxPair({t[OUTL], ExtTauIdx});
      int t2 = G2.AddTidxPair({ExtTauIdx, t[INR]});
      int t3 = G3.AddTidxPair({t[OUTR], ExtTauIdx});
      Gidx.push_back(array<int, 3>({t1, t2, t3}));
    }
  }
};

double sigma2::Evaluate() {
  double Factor = 1.0 / pow(2.0 * PI, D);
  // normalization
  if (Order == 0) {
    return 1.0;
  } else if (Order == 1) {
    double GWeight = Prop.Green(-EPS, Var.LoopMom[0] + Var.LoopMom[1], UP, 0);
    double VerWeight = Prop.Interaction(Var.LoopMom[1], 0);
    // cout << GWeight << ", " << VerWeight << endl;
    return GWeight * VerWeight * Factor;
  }

  // loop order >=2
  vertex4 &Ver4 = Vertex;
  G1.K = Var.LoopMom[1];
  G2.K = Var.LoopMom[2];
  G3.K = Var.LoopMom[0] + Var.LoopMom[2] - Var.LoopMom[1];

  G1.Evaluate();
  G2.Evaluate();
  G3.Evaluate();

  Vertex.Evaluate(Var.LoopMom[0], G1.K, G2.K, G3.K, false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &gidx = Gidx[i];
    double temp = (SPIN * Ver4.Weight[i][DIR] + Ver4.Weight[i][EX]);
    // double temp = Prop.Interaction(G1.K - G2.K, 0);
    // attach four G
    temp *= G1[gidx[1]];
    temp *= G2[gidx[2]];
    temp *= G3[gidx[3]];
    // temp *= Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[0], G1.K, UP, 0);
    // temp *= Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[0], G3.K, UP, 0);
    // temp *= Prop.Green(Var.Tau[0] - Var.Tau[ExtTauIdx], G2.K, UP, 0);
    temp *= Prop.Interaction(G1.K - Var.LoopMom[0], 0);
    if (Order == 3) {
      cout << "(" << Ver4.Tpair[i][0] << ", " << Ver4.Tpair[i][1] << ", "
           << Ver4.Tpair[i][2] << ", " << Ver4.Tpair[i][3]
           << "): " << Ver4.Weight[i][EX] << endl;
      // cout << ExtTauIdx << ", " << Var.Tau[ExtTauIdx] << endl;
    }
    Weight += temp;
  }
  // cout << Order << ", " << Weight << endl;
  return Weight * Factor * Factor;

  // double Factor = 1.0 / pow(2.0 * PI, D);
  // int ExtTauIdx = MaxTauNum - 1;
  // momentum &ExtK = Var.LoopMom[0];
  // momentum &K1 = Var.LoopMom[1];
  // momentum &K2 = Var.LoopMom[2];
  // momentum K3 = K1 + K2 - ExtK;
  // double Tau = Var.Tau[ExtTauIdx] - Var.Tau[0];
  // double G1 = Prop.Green(Tau, K1, UP, 0);
  // double G2 = Prop.Green(Tau, K2, UP, 0);
  // double G3 = Prop.Green(-Tau, K3, UP, 0);
  // double VerDir = Prop.Interaction(K1 - ExtK, 0);
  // double VerEx = -Prop.Interaction(K2 - ExtK, 0);
  // double VerW = 2.0 * VerDir * VerEx;
  // double Weight = VerW * G1 * G2 * G3 * Factor * Factor * 0.5;
  // return Weight;
}