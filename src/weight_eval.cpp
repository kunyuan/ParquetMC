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

double weight::Evaluate(int LoopNum) {
  switch (DiagType) {
  case diagram::GAMMA:
    return EvaluateGamma(LoopNum);
    break;
  case diagram::SIGMA:
    return EvaluateSigma(LoopNum);
    break;
  case diagram::POLAR:
    return EvaluatePolar(LoopNum);
    break;
  default:
    break;
  }
}

double weight::EvaluateGamma(int LoopNum) {
  if (LoopNum == 0) {
    // normalization
    return 1.0;
  } else {
    // if (Channel != dse::T)
    //   return 0.0;
    ver4 &Root = Ver4Root[LoopNum];
    if (Root.Weight.size() != 0) {
      Vertex4(Root, true);
      for (auto &w : ChanWeight)
        // collapse all channel to I
        ChanWeight[0] += w;
      // cout << ChanWeight[0].Abs() << endl;
      return ChanWeight[0].Abs();
    } else
      return 0.0;
  }
}

double weight::EvaluatePolar(int LoopNum) {
  dse::polar &P = Polar[LoopNum];
  ver4 &Root = P.Vertex;
  if (Root.Weight.size() != 0) {
    Vertex4(Root, false);

    // evaluate all possible G
    for (int i = 0; i < 4; ++i)
      EvaluateG(P.G[i], Var.LoopMom[i]);

    int Size = Root.Weight.size();
    P.Weight.SetZero();
    for (int i = 0; i < Size; ++i) {
      auto &Gidx = P.Gidx[i];
      auto &Weight = Root.Weight[i];

      // attach four G
      for (int j = 0; j < 4; ++j)
        Weight *= P.G[j][Gidx[j]].Weight;

      P.Weight += Weight;
    }

    return P.Weight.Abs();
  } else
    return 0.0;
}

double weight::EvaluateSigma(int LoopNum, bool IsFast) {
  dse::sigma &Si = Sigma[LoopNum];
  ver4 &Root = Si.Vertex;
  if (Root.Weight.size() != 0) {
    Vertex4(Root, false);
    EvaluateG(Si.G, *Si.K);
    int Size = Root.Weight.size();
    for (int i = 0; i < Size; ++i) {
      int Gidx = Si.Gidx[i];
      auto &Weight = Root.Weight[i];
      if (IsFast)
        Si.Weight[0] += Weight.Sum() * Si.G[Gidx].Weight;
      else
        Si.Weight[Si.T[i]] += Weight.Sum() * Si.G[Gidx].Weight;
    }
    return Si.Weight[0];
  } else
    return 0.0;
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
      for (auto &w : ChanWeight)
        w.SetZero();
    }

    ChanUST(Ver4, IsFast);
    if (Ver4.LoopNum >= 3)
      ChanI(Ver4, IsFast);
  }
  return;
}

void weight::EvaluateG(vector<green> &G, const momentum &K) {
  // calculate G table
  for (auto &g : G) {
    g.Weight = Fermi.Green(Var.Tau[g.T[OUT]] - Var.Tau[g.T[IN]], K, UP, 0, 0);
    // cout << g.Weight << endl;
  }
}

void weight::ChanUST(ver4 &Ver4, bool IsFast) {
  double Weight = 0.0;
  double Ratio, dTau, ProjFactor;
  double DirW, ExW;
  const auto &LegK = Ver4.LegK;
  auto &G = Ver4.G;
  double Factor = 1.0 / pow(2.0 * PI, D);

  // if (Ver4.Level == 0 && Ver4.LoopNum == 2)
  //   cout << "\nNew Ver4" << endl;

  // calculate K table
  Ver4.K[0] = Var.LoopMom[Ver4.Loopidx];
  EvaluateG(G[0], Ver4.K[0]);
  for (auto chan : Ver4.Channel) {
    // cout << "chan: " << chan << endl;
    switch (chan) {
    case T:
      Ver4.K[T] = *LegK[OUTL] + Ver4.K[0] - *LegK[INL];
      if (Ver4.InBox)
        EvaluateG(G[chan], Ver4.K[0]);
      else
        EvaluateG(G[chan], Ver4.K[chan]);
      break;
    case U:
      // ASSERT_ALLWAYS(Ver4.InBox == false, "fail test");
      Ver4.K[U] = *LegK[OUTR] + Ver4.K[0] - *LegK[INL];
      EvaluateG(G[chan], Ver4.K[chan]);
      break;
    case S:
      // ASSERT_ALLWAYS(Ver4.InBox == false, "fail test");
      Ver4.K[S] = *LegK[INL] + *LegK[INR] - Ver4.K[0];
      EvaluateG(G[chan], Ver4.K[chan]);
      break;
    case TC:
      // Ver4.K[T] may have not yet been initialized before this step
      Ver4.K[T] = *LegK[OUTL] + Ver4.K[0] - *LegK[INL];
      EvaluateG(G[chan], Ver4.K[0]);
      break;
    case UC:
      // Ver4.K[U] may have not yet been initialized before this step
      Ver4.K[U] = *LegK[OUTR] + Ver4.K[0] - *LegK[INL];
      EvaluateG(G[chan], Ver4.K[0]);
      break;
    }
  }
  // cout << "InL=(" << (*LegK[INL])[0] << ", " << (*LegK[INL])[1] << ", "
  //      << (*LegK[INL])[2] << ")" << endl;
  // cout << "OutL=(" << (*LegK[OUTL])[0] << ", " << (*LegK[OUTL])[1] << ", "
  //      << (*LegK[OUTL])[2] << ")" << endl;
  // cout << "InR=(" << (*LegK[INR])[0] << ", " << (*LegK[INR])[1] << ", "
  //      << (*LegK[INR])[2] << ")" << endl;
  // cout << "OutR=(" << (*LegK[OUTR])[0] << ", " << (*LegK[OUTR])[1] << ", "
  //      << (*LegK[OUTR])[2] << ")" << endl;
  // cout << "KU=(" << (Ver4.K[U])[0] << ", " << (Ver4.K[U])[1] << ", "
  //      << (Ver4.K[U])[2] << ")" << endl;
  // cout << Ver4.K[UC].norm() << ", " << Ver4.K[0].norm() << endl;
  // cout << "LegK0: " << (*Ver4.LegK[0]).norm() << endl;
  // cout << "LegK1: " << (*Ver4.LegK[1]).norm() << endl;
  // cout << "LegK2: " << (*Ver4.LegK[2]).norm() << endl;
  // cout << "LegK3: " << (*Ver4.LegK[3]).norm() << endl;
  // cout << Ver4.LegK[0] << ", " << Ver4.LegK[1] << ", " << Ver4.LegK[2] << ",
  // "
  //      << Ver4.LegK[3] << endl;
  // cout << "K0: " << Ver4.K[0].norm() << ", " << &Ver4.K[0] << endl;

  ///////////// Check if the projected counter-terms exist or not ///////

  // for vertex4 with one or more loops
  for (auto &b : Ver4.Bubble) {
    ver4 &LVer = b.LVer;
    ver4 &RVer = b.RVer;
    Vertex4(LVer, IsFast);
    Vertex4(RVer, IsFast);

    ProjFactor = SymFactor[b.Channel] * Factor;
    if (b.Channel == TC || b.Channel == UC || Ver4.InBox)
      ProjFactor *= Para.Lambda / (8.0 * PI * Para.Nf);
    // else
    //   ProjFactor = 0.0;

    for (auto &map : b.Map) {
      Weight =
          ProjFactor * G[0][map.G0idx].Weight * G[b.Channel][map.Gidx].Weight;

      // if (b.Channel == UC || b.Channel == TC || Ver4.InBox) {
      //   ASSERT_ALLWAYS(
      //       abs(G[0][map.G0idx].Weight * G[b.Channel][map.Gidx].Weight -
      //           Fermi.Green(Para.Beta / 2.0, Ver4.K[0], UP, 0, 0) *
      //               Fermi.Green(-Para.Beta / 2.0, Ver4.K[0], UP, 0, 0)) <
      //           1.0e-10,
      //       "GG does not match!");
      // }

      auto &Lw = LVer.Weight[map.LVerTidx];
      auto &Rw = RVer.Weight[map.RVerTidx];

      // if (b.Channel == UC || b.Channel == TC || Ver4.InBox) {
      //   ASSERT_ALLWAYS(abs(Lw[EX]) < 1.0e-10, "LWrong! " << Lw[EX]);
      //   ASSERT_ALLWAYS(abs(Rw[EX]) < 1.0e-10, "RWrong! " << Rw[EX]);
      // }

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