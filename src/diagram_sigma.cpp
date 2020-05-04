#include "diagram.h"

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void sigma::Build(int order) {
  ASSERT_ALLWAYS(order >= 1, "Sigma LoopNum must be larger than 0!");
  Order = order;

  Ver4.Level = 0;
  Ver4.LoopNum = LoopNum - 1;
  Ver4.TauNum = LoopNum;
  Ver4.Side = RIGHT;
  Ver4.InTL = 0;
  Ver4.Loopidx = 2;
  Ver4.InBox = false;
  Ver4.Channel = {dse::T, dse::TC};
  int InTL = Ver4.InTL;

  vector<channel> FULL = {I, T, U, S, TC, UC};
  vector<channel> F = {TC, UC}; // the bare interaction is automatically
                                // included
  vector<channel> FULL_CT = {I, T, TC};
  vector<channel> F_CT = {TC}; // the bare interaction is automatically included

  for (auto &chan : Ver4.Channel) {

    for (int ol = 0; ol < Ver4.LoopNum; ol++) {
      bubble bub;
      ////////////////////   Right SubVer  ///////////////////
      int Level = Ver4.Level + 1;
      int oR = Ver4.LoopNum - 1 - ol;
      int RInTL = Ver4.InTL + (ol + 1);
      int Llopidx = Ver4.Loopidx + 1;
      int Rlopidx = Ver4.Loopidx + 1 + ol;

      switch (chan) {
      case T:
        bub.LVer = Factory.Vertex(Level, ol, Llopidx, InTL, F, LEFT, false);
        bub.RVer =
            Factory.Vertex(Level, oR, Rlopidx, RInTL, FULL, RIGHT, false);
        break;
      case TC:
        // continue;
        bub.LVer = Factory.Vertex(Level, ol, Llopidx, InTL, F_CT, LEFT, true);
        bub.RVer =
            Factory.Vertex(Level, oR, Rlopidx, RInTL, FULL_CT, RIGHT, true);
        break;
      default:
        ABORT("The channel does not exist!");
        break;
      }

      bub.Channel = chan;
      bub.Map = CreateIndexMap(Ver4, bub.LVer, bub.RVer, chan);
      Ver4.Bubble.push_back(bub);
    }
  }

  // create the G List
  for (int i = 0; i < Ver4.T.size(); ++i) {
    auto &t = Ver4.T[i];
    // the G is from OUTL to INR
    int Gidx = AddToGList(Sigma.G, {t[OUTL], t[INR]});
    // the tau of Sigma is from 0 to OUTR
    int SigmaTidx = AddToSingleTList(Sigma.T, t[OUTR]);

    Sigma.Gidx.push_back(Gidx);
    Sigma.SigTidx.push_back(SigmaTidx);
  }

  Ver4.Weight.resize(Ver4.T.size());
  Sigma.Weight.resize(Sigma.T.size());
  return Sigma;
}