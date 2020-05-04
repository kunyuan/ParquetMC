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
    Bubble.push_back(bubble());
    auto &bub = Bubble.back();

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
    }
  }
}