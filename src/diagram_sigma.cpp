#include "diagram.h"

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void sigma::Build(int order) {
  ASSERT_ALLWAYS(order >= 0, "Sigma LoopNum must be larger than 0!");
  Order = order;
  if (order <= 1)
    return;

  vector<channel> FULL = {I, T, U, S, TC, UC};
  vector<channel> F = {TC, UC};
  // the bare interaction is automatically included

  for (int ol = 0; ol < LoopNum() - 1; ol++) {
    bubble bub;

    ////////////////////   Right SubVer  ///////////////////
    int lvl = 0;
    int oR = LoopNum() - 2 - ol;
    int LInTL = 0;
    int RInTL = LInTL + (ol + 1);
    int Llopidx = 3; // ExtK: 0, G1: 1, G2: 2
    int Rlopidx = 3 + ol;

    ASSERT(RInTL == TauNum() - 1, "RInTL must be the last tau idx!");

    bub.LVer.Build(lvl, ol, Llopidx, LInTL, FULL, LEFT, false);
    bub.RVer.Build(lvl, oR, Rlopidx, RInTL, F, RIGHT, false);

    for (int ol = 0; ol < bub.LVer.Tpair.size(); ++ol) {
      // Assume the RVer is equal-time!
      auto &t = bub.LVer.Tpair[ol];

      int G1idx = G1.AddTidxPair({t[OUTL], RInTL});
      int G2idx = G2.AddTidxPair({t[OUTR], RInTL});
      int G3idx = G3.AddTidxPair({RInTL, t[INR]});

      Map.push_back({ol, G1idx, G2idx, G3idx});
    }

    Bubble.push_back(bub);
  }
}

double sigma::Evaluate() {
  // normalization
  double Factor = 1.0 / pow(2.0 * PI, D);
  if (Order == 0) {
    return 1.0;
  } else if (Order == 1) {
    double GWeight = Prop.Green(-EPS, Var.LoopMom[0] + Var.LoopMom[1], UP, 0);
    double VerWeight = Prop.Interaction(Var.LoopMom[1], 0);
    // cout << GWeight << ", " << VerWeight << endl;
    return GWeight * VerWeight * Factor;
  }

  // Sigma with LoopNum>=2
  G1.K = Var.LoopMom[1];
  G2.K = Var.LoopMom[2];
  G3.K = Var.LoopMom[1] + Var.LoopMom[2] - Var.LoopMom[0];
  G1.Evaluate();
  G2.Evaluate();
  G3.Evaluate();

  double Weight = 0.0;

  for (auto &b : Bubble) {
    b.LVer.Evaluate(Var.LoopMom[0], G1.K, G3.K, G2.K, false);
    b.RVer.Evaluate(G2.K, G3.K, G1.K, Var.LoopMom[0], false);

    for (auto &map : Map) {
      auto &LvW = b.LVer.Weight[map[0]];
      auto &RvW = b.RVer.Weight[map[0]];
      double w = (LvW[DIR] * RvW[DIR] + LvW[EX] * RvW[EX]) * SPIN;
      w += LvW[DIR] * RvW[EX] + LvW[EX] + RvW[DIR];
      w *= G1[map[1]] * G2[map[2]] * G3[map[3]] * 0.5;
      Weight += w;
    }
  }
  return Weight * Factor * Factor;
}