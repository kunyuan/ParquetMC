#include "diagram.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <string>

using namespace dse;
using namespace std;

dse::polar dse::BuildPolar(int LoopNum, array<momentum *, 4> ExtLegK) {
  dse::polar Polar;
  Polar.LoopNum = LoopNum;
  Polar.TauNum = LoopNum + 1;

  vector<channel> Chan = {I, T, U, S, TC, UC};
  // vector<channel> Chan = {T, TC};
  // vector<channel> Chan = {U, UC};
  verDiag Factory;
  Polar.Vertex = Factory.Vertex(0,           // level
                                LoopNum - 2, // loopNum
                                4, // loop index of the first internal K
                                2, // tau index of the InTL leg
                                Chan, RIGHT, false);
  Factory.ResetMomMap(Polar.Vertex, ExtLegK);
  for (auto &t : Polar.Vertex.T) {
    AddToGList(Polar.G[INL], {0, t[INL]});
    AddToGList(Polar.G[OUTL], {t[OUTL], 0});
    AddToGList(Polar.G[INR], {1, t[INR]});
    AddToGList(Polar.G[OUTR], {t[OUTR], 1});
  }
  return Polar;
};

dse::sigma dse::BuildSigma(int LoopNum, momentum *ExtK) {
  dse::sigma Sigma;
  Sigma.LoopNum = LoopNum;
  Sigma.TauNum = LoopNum;

  auto &Ver4 = Sigma.Vertex;

  verDiag Factory;

  Ver4.Level = 0;
  Ver4.LoopNum = LoopNum - 1;
  Ver4.TauNum = LoopNum;
  Ver4.Side = RIGHT;
  Ver4.InTL = 0;
  Ver4.Loopidx = 0;
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
  /////////////////// set K table ////////////////////
  auto &K = Ver4.K;
  auto &LegK = Ver4.LegK;
  Ver4.LegK[INL] = ExtK;
  Ver4.LegK[OUTR] = ExtK;
  Ver4.LegK[OUTL] = &Sigma.K;
  Ver4.LegK[INR] = &Sigma.K;
  array<momentum *, 4> LLegK, RLegK;

  for (auto &bub : Ver4.Bubble) {
    auto chan = bub.Channel;
    ASSERT_ALLWAYS(chan == T || chan == TC, "Should be T or TC channel!");

    LLegK = {LegK[INL], LegK[OUTL], &K[T], &K[0]};
    RLegK = {&K[0], &K[T], LegK[INR], LegK[OUTR]};
    Factory.ResetMomMap(bub.LVer, LLegK);
    Factory.ResetMomMap(bub.RVer, RLegK);
  }

  // create the G List
  for (auto &t : Ver4.T) {
    // the G is from OUTL to INR
    AddToGList(Sigma.G, {t[OUTL], t[INR]});
    // the tau of Sigma is from 0 to OUTR
    Sigma.T.push_back(t[OUTR]);
  }

  // for (auto &g : Sigma.G)
  //   cout << g.T[IN] << "->" << g.T[OUT] << endl;

  Sigma.Weight.resize(Sigma.T.size());
  return Sigma;
}