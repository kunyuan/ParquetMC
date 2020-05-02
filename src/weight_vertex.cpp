#include "diagram.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include "utility/vector.h"
#include "weight.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

void weight::Vertex4(ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
                     const momentum &KInR, const momentum &KOutR, bool IsFast) {
  if (Ver4.LoopNum == 0) {
    Ver0(Ver4, KInL, KOutL, KInR, KOutR);
  } else {
    if (Ver4.Level > 0 || IsFast == false) {
      for (auto &w : Ver4.Weight)
        w.SetZero();
    } else {
      for (auto &w : ChanWeight)
        w.SetZero();
    }

    ChanUST(Ver4, KInL, KOutL, KInR, KOutR, IsFast);
    if (Ver4.LoopNum >= 3)
      ChanI(Ver4, KInL, KOutL, KInR, KOutR, IsFast);
  }
  return;
}

void weight::Ver0(ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
                  const momentum &KInR, const momentum &KOutR) {
  // only bare coupling
  if (DiagType == POLAR)
    Ver4.Weight[0] = VerQTheta.Interaction(KInL, KOutL, KInR, KOutR, 0.0, false,
                                           Ver4.InBox, Var.LoopMom[0].norm());
  else
    Ver4.Weight[0] =
        VerQTheta.Interaction(KInL, KOutL, KInR, KOutR, 0.0, false, Ver4.InBox);
  // if (abs(Ver4.Weight[0][DIR] + 8.377) > 0.1) {
  //   cout << "Ver0 " << Ver4.Weight[0][DIR] << endl;
  //   ABORT("Err");
  // }
  // if (Ver4.Level == 1) {
  //   cout << "ver0: " << Ver4.Weight[0][DIR] << Ver4.Weight[0][EX] << endl;
  //   cout << "ver0 LegK0: " << (*Ver4.LegK[0]).norm() << endl;
  //   cout << "ver0 LegK1: " << (*Ver4.LegK[1]).norm() << endl;
  //   cout << "ver0 LegK2: " << (*Ver4.LegK[2]).norm() << endl;
  //   cout << "ver0 LegK3: " << (*Ver4.LegK[3]).norm() << endl;
  //   cout << Ver4.LegK[0] << ", " << Ver4.LegK[1] << ", " << Ver4.LegK[2] <<
  //   ", "
  //        << Ver4.LegK[3] << endl;
  // }
  return;
}

void weight::EvaluateG(vector<green> &G, const momentum &K) {
  // calculate G table
  for (auto &g : G) {
    g.Weight = Fermi.Green(Var.Tau[g.T[OUT]] - Var.Tau[g.T[IN]], K, UP, 0, 0);
    // cout << g.Weight << endl;
  }
}

void weight::ChanUST(ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
                     const momentum &KInR, const momentum &KOutR, bool IsFast) {
  double Weight = 0.0;
  double Ratio, dTau, ProjFactor;
  double DirW, ExW;
  // const auto &LegK = Ver4.LegK;
  auto &G = Ver4.G;
  auto &K = Ver4.K;
  double Factor = 1.0 / pow(2.0 * PI, D);

  // calculate K table
  auto &K0 = Var.LoopMom[Ver4.Loopidx];
  EvaluateG(G[0], K0);
  for (auto chan : Ver4.Channel) {
    switch (chan) {
    case T:
      K[T] = KOutL + K0 - KInL;
      if (Ver4.InBox)
        EvaluateG(G[chan], K0);
      else
        EvaluateG(G[chan], Ver4.K[chan]);
      break;
    case U:
      K[U] = KOutR + K0 - KInL;
      EvaluateG(G[chan], Ver4.K[chan]);
      break;
    case S:
      K[S] = KInL + KInR - K0;
      EvaluateG(G[chan], Ver4.K[chan]);
      break;
    case TC:
      // Ver4.K[T] may have not yet been initialized before this step
      K[T] = KOutL + K0 - KInL;
      EvaluateG(G[chan], K0);
      break;
    case UC:
      // Ver4.K[U] may have not yet been initialized before this step
      K[U] = KOutR + K0 - KInL;
      EvaluateG(G[chan], K0);
      break;
    }
  }
  ///////////// Check if the projected counter-terms exist or not ///////
  // for vertex4 with one or more loops
  for (auto &b : Ver4.Bubble) {
    ver4 &LVer = b.LVer;
    ver4 &RVer = b.RVer;

    if (b.Channel == T || b.Channel == TC) {
      Vertex4(LVer, KInL, KOutL, K[T], K0, IsFast);
      Vertex4(RVer, K0, K[T], KInR, KOutR, IsFast);
    } else if (b.Channel == U || b.Channel == UC) {
      Vertex4(LVer, KInL, KOutR, K[U], K0, IsFast);
      Vertex4(RVer, K0, K[U], KInR, KOutL, IsFast);
    } else if (b.Channel == S) {
      Vertex4(LVer, KInL, K[S], KInR, K0, IsFast);
      Vertex4(RVer, K0, KOutL, K[S], KOutR, IsFast);
    }

    ProjFactor = SymFactor[b.Channel] * Factor;
    if (b.Channel == TC || b.Channel == UC || Ver4.InBox)
      ProjFactor *= Para.Lambda / (8.0 * PI * Para.Nf);
    // else
    //   ProjFactor = 0.0;

    for (auto &map : b.Map) {
      Weight =
          ProjFactor * G[0][map.G0idx].Weight * G[b.Channel][map.Gidx].Weight;

      auto &Lw = LVer.Weight[map.LVerTidx];
      auto &Rw = RVer.Weight[map.RVerTidx];

      if (b.Channel == T || b.Channel == TC) {
        DirW = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        ExW = Lw[EX] * Rw[EX];
      } else if (b.Channel == U || b.Channel == UC) {
        ExW = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        DirW = Lw[EX] * Rw[EX];
      } else if (b.Channel == S) {
        // see the note "code convention"
        DirW = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        ExW = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX];
      }

      if (Ver4.Level > 0 || IsFast == false) {
        Ver4.Weight[map.Tidx][DIR] += DirW * Weight;
        Ver4.Weight[map.Tidx][EX] += ExW * Weight;
      } else {
        channel chan = ChanMap[b.Channel];
        ChanWeight[chan][DIR] += DirW * Weight;
        ChanWeight[chan][EX] += ExW * Weight;
        // cout << DirW << ", " << ExW << ", " << Weight << ", " << ProjFactor
        //      << endl;
      }
    }
  }
}

void weight::ChanI(dse::ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
                   const momentum &KInR, const momentum &KOutR, bool IsFast) {
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