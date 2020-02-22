#include "weight.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

void weight::Initialization() {

  array<momentum *, 4> ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1],
                                  &Var.LoopMom[2], &Var.LoopMom[3]};

  vector<channel> Chan = {I, T, U, S, TC, UC};
  for (int order = 1; order <= Para.Order; order++) {
    LOG_INFO("Generating order " << order);
    Ver4Root[order] = VerDiag.Vertex(0,     // level
                                     order, // loopNum
                                     4, // loop index of the first internal K
                                     0, // tau index of the InTL leg
                                     Chan, ExtLegK, RIGHT, false);
    if (order < 5)
      LOG_INFO(VerDiag.ToString(Ver4Root[order]));
  }
  for (int order = 1; order <= Para.Order; order++) {
    LOG_INFO(fmt::format(" Order {0}: ExtT Num {1}\n", order,
                         Ver4Root[order].T.size()));
  }
}

double weight::Evaluate(int LoopNum, diagram Diagram) {
  if (LoopNum == 0) {
    // normalization
    return 1.0;
  } else {
    // if (Channel != dse::T)
    //   return 0.0;
    ver4 &Root = Ver4Root[LoopNum];
    if (Root.Weight.size() != 0) {
      Vertex4(Root, true);
      for (auto &w : _ChanWeight)
        _ChanWeight[0] += w;
      return _ChanWeight[0].Abs();
    }
  }
}

void weight::Ver0(ver4 &Ver4) {
  // only bare coupling
  Ver4.Weight[0] = VerQTheta.Interaction(Ver4.LegK, 0.0, false, Ver4.InBox);
  return;
}
void weight::Vertex4(ver4 &Ver4, bool IsFast) {
  // cout << Ver4.LoopNum << endl;
  if (Ver4.LoopNum == 0) {
    Ver0(Ver4);
  } else {
    if (Ver4.Level > 0 || IsFast == false) {
      for (auto &w : Ver4.Weight)
        w.SetZero();
    } else {
      for (auto &w : _ChanWeight)
        w.SetZero();
    }

    ChanUST(Ver4, IsFast);
    if (Ver4.LoopNum >= 3)
      ChanI(Ver4, IsFast);
  }
  return;
}

void weight::ChanUST(ver4 &Ver4, bool IsFast) {
  double Weight = 0.0;
  double Ratio, dTau, ProjFactor;
  double DirW, ExW;
  const auto &LegK = Ver4.LegK;
  auto &G = Ver4.G;
  double Factor = 1.0 / pow(2.0 * PI, D * Ver4.LoopNum);

  // calculate K table
  Ver4.K[0] = Var.LoopMom[Ver4.Loopidx];
  for (auto chan : Ver4.Channel) {
    switch (chan) {
    case T:
      Ver4.K[T] = *LegK[OUTL] + Ver4.K[0] - *LegK[INL];
      break;
    case U:
      Ver4.K[U] = *LegK[OUTR] + Ver4.K[0] - *LegK[INL];
      break;
    case S:
      Ver4.K[S] = *LegK[INL] + *LegK[INR] - Ver4.K[0];
      break;
    case TC:
    case UC:
      ProjFactor = Para.Lambda / (8.0 * PI * Para.Nf);
      Ver4.K[chan] = Ver4.K[0];
      break;
    }

    // calculate G table
    for (auto &g : G[chan]) {
      dTau = Var.Tau[g.T[OUT]] - Var.Tau[g.T[IN]];
      g.Weight = Fermi.Green(dTau, Ver4.K[chan], UP, 0, 0);
    }
  }
  ///////////// Check if the projected counter-terms exist or not ///////

  // for vertex4 with one or more loops
  for (auto &b : Ver4.Bubble) {
    ver4 &LVer = b.LVer;
    ver4 &RVer = b.RVer;
    Vertex4(LVer, IsFast);
    Vertex4(RVer, IsFast);

    ProjFactor = SymFactor[b.Channel] * Factor;
    if (b.Channel == TC || b.Channel == UC)
      ProjFactor *= Para.Lambda / (8.0 * PI * Para.Nf);

    for (auto &map : b.Map) {
      Weight = ProjFactor;
      Weight *= G[0][map.G0idx].Weight * G[b.Channel][map.Gidx].Weight;

      auto &Lw = LVer.Weight[map.LVerTidx];
      auto &Rw = RVer.Weight[map.RVerTidx];

      switch (b.Channel) {
      case T:
      case TC:
        DirW = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        ExW = Lw[EX] * Rw[EX];
        break;
      case U:
      case UC:
        ExW = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        DirW = Lw[EX] * Rw[EX];
        break;
      case S:
        // see the note "code convention"
        DirW = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        ExW = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX];
        break;
      }

      if (Ver4.Level > 0) {
        Ver4.Weight[map.Tidx][DIR] += DirW * Weight;
        Ver4.Weight[map.Tidx][EX] += ExW * Weight;
      } else {
        if (IsFast) {
          // calculate contributions from different channels for the root
          // vertex4
          channel chan = ChanMap[b.Channel];
          _ChanWeight[chan][DIR] += DirW * Weight;
          _ChanWeight[chan][EX] += ExW * Weight;
        } else {
          // TODO: add slow measure
        }
      }
    }
  }
}

void weight::ChanI(dse::ver4 &Ver4, bool IsFast) {
  if (Ver4.LoopNum != 3)
    return;
  // for (auto &Env : Ver4.Envelope) {
  //   const momentum &InL = *Env.LegK[INL];
  //   const momentum &OutL = *Env.LegK[OUTL];
  //   const momentum &InR = *Env.LegK[INR];
  //   const momentum &OutR = *Env.LegK[OUTR];

  //   auto &G = Env.G;

  //   *G[3].K = *G[0].K + *G[1].K - InL;
  //   *G[4].K = *G[1].K + *G[2].K - OutL;
  //   *G[5].K = *G[0].K + InR - *G[2].K;
  //   *G[6].K = *G[1].K + *G[2].K - OutR;
  //   *G[7].K = *G[2].K + OutR - *G[0].K;
  //   *G[8].K = *G[2].K + OutL - *G[0].K;

  //   for (auto &g : Env.G)
  //     g.Weight = Fermi.Green(Var.Tau[g.OutT] - Var.Tau[g.InT], *(g.K), UP,
  //     0,
  //                            Var.CurrScale);

  //   for (auto &subVer : Env.Ver)
  //     Vertex4(subVer);

  //   double Weight = 0.0;
  //   double ComWeight = 0.0;
  //   for (auto &map : Env.Map) {
  //     auto &SubVer = Env.Ver;
  //     auto &GT = map.GT;
  //     auto &G = Env.G;
  //     ComWeight = G[0].Weight * G[1].Weight * G[2].Weight * G[3].Weight;
  //     // cout << "G: " << ComWeight << endl;
  //     ComWeight *= SubVer[0].Weight[map.LDVerTidx];
  //     // cout << "Ver: " << SubVer[0].Weight[map.LDVerT] << endl;
  //     // cout << "T: " << map.LDVerT << endl;

  //     Weight = Env.SymFactor[0] * ComWeight;
  //     Weight *= SubVer[1].Weight[map.LUVerTidx];
  //     Weight *= SubVer[3].Weight[map.RDVerTidx];
  //     Weight *= SubVer[6].Weight[map.RUVerTidx];
  //     Weight *= G[4].Weight * G[5].Weight;
  //     Ver4.Weight[map.Tidx[0]] += Weight;

  //     Weight = Env.SymFactor[1] * ComWeight;
  //     Weight *= SubVer[2].Weight[map.LUVerTidx];
  //     Weight *= SubVer[3].Weight[map.RDVerTidx];
  //     Weight *= SubVer[7].Weight[map.RUVerTidx];
  //     Weight *= G[6].Weight * G[5].Weight;
  //     Ver4.Weight[map.Tidx[1]] += Weight;
  //     // cout << Weight << endl;

  //     Weight = Env.SymFactor[2] * ComWeight;
  //     Weight *= SubVer[1].Weight[map.LUVerTidx];
  //     Weight *= SubVer[4].Weight[map.RDVerTidx];
  //     Weight *= SubVer[8].Weight[map.RUVerTidx];
  //     Weight *= G[4].Weight * G[7].Weight;
  //     Ver4.Weight[map.Tidx[2]] += Weight;
  //     // cout << Weight << endl;

  //     Weight = Env.SymFactor[3] * ComWeight;
  //     Weight *= SubVer[2].Weight[map.LUVerTidx];
  //     Weight *= SubVer[5].Weight[map.RDVerTidx];
  //     Weight *= SubVer[9].Weight[map.RUVerTidx];
  //     Weight *= G[6].Weight * G[8].Weight;
  //     Ver4.Weight[map.Tidx[3]] += Weight;
  //     // cout << Weight << endl;

  //     // if (map.LDVerT == 0 && map.LUVerT == 0 && map.RDVerT == 0 &&
  //     //     map.RUVerT == 0) {
  //     // cout << "Com: " << ComWeight << endl;
  //     // cout << "G[4]: " << G[4](GT[4]) << endl;
  //     // cout << "G[5]: " << G[5](GT[5]) << endl;
  //     // cout << SubVer[1].Weight[map.LUVerT] << endl;
  //     // cout << SubVer[3].Weight[map.RDVerT] << endl;
  //     // cout << SubVer[6].Weight[map.RUVerT] << endl;
  //     // cout << "First: " << Weight << endl;
  //     // }
  //   }
  // }

  return;
}
