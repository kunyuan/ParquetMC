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
    // vector<channel> Chan = {I, T, U, S, TC, UC};
    // vector<channel> Chan = {T, U, S, I};
    vector<channel> Chan = {T};
    // vector<channel> Chan = {
    //     T,
    // };
    Vertex.Build(0,         // level
                 Order - 2, // loopNum
                 3,         // loop index of the first internal K of the vertex
                 0,         // tau index of the InTL leg
                 Chan, RIGHT);

    for (auto &t : Vertex.Tpair) {
      // cout << t[INL] << ": " << t[OUTL] << ": " << t[INR] << ": " << t[OUTR]
      //      << endl;
      int t1 = G1.AddTidxPair({t[OUTL], ExtTauIdx});
      int t2 = G2.AddTidxPair({ExtTauIdx, t[INR]});
      int t3 = G3.AddTidxPair({t[OUTR], ExtTauIdx});
      // cout << t1 << "; " << t2 << "; " << t3 << endl;
      Gidx.push_back(array<int, 3>({t1, t2, t3}));
    }
    // cout << G1._Tpair[0][0] << G1._Tpair[0][1] << endl;
    // cout << G2._Tpair[0][0] << G2._Tpair[0][1] << endl;
    // cout << G3._Tpair[0][0] << G3._Tpair[0][1] << endl;
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
  } // loop order >=2
  vertex4 &Ver4 = Vertex;
  G1.K = Var.LoopMom[1];
  G2.K = Var.LoopMom[2];
  G3.K = Var.LoopMom[0] + Var.LoopMom[2] - Var.LoopMom[1];

  // cout << "Evaluate G1, 2, 3" << endl;
  G1.Evaluate();
  G2.Evaluate();
  G3.Evaluate();
  // cout << G1._Tpair[0][0] << G1._Tpair[0][1] << endl;
  // cout << G2._Tpair[0][0] << G2._Tpair[0][1] << endl;
  // cout << G3._Tpair[0][0] << G3._Tpair[0][1] << endl;

  Vertex.Evaluate(Var.LoopMom[0], G1.K, G2.K, G3.K, false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &gidx = Gidx[i];
    double temp = SPIN * Ver4.Weight[i][DIR] + Ver4.Weight[i][EX];
    temp *= G1[gidx[0]];
    temp *= G2[gidx[1]];
    temp *= G3[gidx[2]];
    temp *= Prop.Interaction(G1.K - Var.LoopMom[0], 0);
    Weight += temp;
  }
  // cout << Order << ", " << Weight << endl;
  return Weight * Factor * Factor;
}

bool sigma2::Test() {
  return false;
  // Two-loop sigma
  if (Order != 3)
    return false;

  double Factor = 1.0 / pow(2.0 * PI, D);
  int ExtTauIdx = MaxTauNum - 1;
  momentum &K0 = Var.LoopMom[0];
  momentum &K1 = Var.LoopMom[1];
  momentum &K2 = Var.LoopMom[2];
  momentum K3 = K0 + K2 - K1;
  momentum K4 = Var.LoopMom[3];
  double Ver0 = Prop.Interaction(K0 - K1, 0);
  double Common = Ver0 * Factor * Factor * Factor;

  //////    T diagram  ///////////////////

  momentum K5 = K1 + K4 - K0;
  double g1 = Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[0], K1, UP, 0);
  double g2 = Prop.Green(Var.Tau[1] - Var.Tau[ExtTauIdx], K2, UP, 0);
  double g3 = Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[1], K3, UP, 0);
  double g4 = Prop.Green(Var.Tau[1] - Var.Tau[0], K4, UP, 0);
  double g5 = Prop.Green(Var.Tau[0] - Var.Tau[1], K5, UP, 0);
  double Ver1 = Prop.Interaction(K0 - K4, 0);
  double Ver2 = Prop.Interaction(K3 - K4, 0);
  double WeightT = -Ver1 * Ver2 * g1 * g2 * g3 * g4 * g5 * Common;

  // cout << "Test G: " << g1 << ", " << g2 << ", " << g3 << ", " << g4 << ", "
  //      << g5 << endl;
  // cout << "Test Ver: " << Ver1 << ", " << Ver2 << endl;

  ////////   U diagram //////////////////////
  K5 = K4 + K3 - K0;
  g1 = Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[1], K1, UP, 0);
  g2 = Prop.Green(Var.Tau[1] - Var.Tau[ExtTauIdx], K2, UP, 0);
  g3 = Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[0], K3, UP, 0);
  g4 = Prop.Green(Var.Tau[1] - Var.Tau[0], K4, UP, 0);
  g5 = Prop.Green(Var.Tau[0] - Var.Tau[1], K5, UP, 0);
  double Ver1d = Prop.Interaction(K0 - K3, 0);
  double Ver1e = -Prop.Interaction(K0 - K4, 0);
  double Ver2d = Prop.Interaction(K2 - K1, 0);
  double Ver2e = -Prop.Interaction(K2 - K5, 0);
  double WeightU = (SPIN * Ver1d * Ver2d + Ver1d * Ver2e + Ver1e * Ver2d +
                    Ver1e * Ver2e * SPIN) *
                   g1 * g2 * g3 * g4 * g5 * Common;

  ////////   S diagram //////////////////////
  K5 = K2 + K0 - K4;
  g1 = Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[1], K1, UP, 0);
  g2 = Prop.Green(Var.Tau[0] - Var.Tau[ExtTauIdx], K2, UP, 0);
  g3 = Prop.Green(Var.Tau[ExtTauIdx] - Var.Tau[1], K3, UP, 0);
  g4 = Prop.Green(Var.Tau[1] - Var.Tau[0], K4, UP, 0);
  g5 = Prop.Green(Var.Tau[1] - Var.Tau[0], K5, UP, 0);
  Ver1d = Prop.Interaction(K0 - K5, 0);
  Ver1e = -Prop.Interaction(K0 - K4, 0);
  Ver2d = Prop.Interaction(K1 - K4, 0);
  Ver2e = -Prop.Interaction(K1 - K5, 0);
  double WeightS = -(Ver1d * Ver2d + SPIN * Ver1d * Ver2e +
                     SPIN * Ver1e * Ver2d + Ver1e * Ver2e) *
                   g1 * g2 * g3 * g4 * g5 * Common * 0.5;

  // cout << "Test G: " << g1 << ", " << g2 << ", " << g3 << ", " << g4 << ", "
  //      << g5 << endl;
  // cout << "Test Ver1: " << Ver1d << ", " << Ver1e << endl;
  // cout << "Test Ver2: " << Ver2d << ", " << Ver2e << endl;
  // cout << "Test Total Ver: "
  //      << (Ver1d * Ver2d + SPIN * Ver1d * Ver2e + SPIN * Ver1e * Ver2d +
  //          Ver1e * Ver2e) *
  //             0.5 * Factor * g4 * g5
  //      << endl;
  double Weight = WeightT + WeightU + WeightS;
  // double Weight = WeightT;
  ASSERT_ALLWAYS(IsEqual(Weight, Evaluate()),
                 "Sigma weight error: " << Weight << " vs " << Evaluate());
  // cout << "G_test=" << G1 << ", " << G2 << ", " << G3 << endl;
  // cout << "Ver_test=" << VerWeightDir << ", " << VerWeightExLeft << endl;

  //   cout << "DIR=" << Weight1 + 2 * Weight2 << ", EX=" << -Weight3 << endl;
  //   return Weight1 + 2 * Weight2;
  // cout << "Pass" << endl;
  return true;
}