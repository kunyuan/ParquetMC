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
      Vertex4(Root,
              Var.LoopMom[INL],  // KInL
              Var.LoopMom[OUTL], // KOutL
              Var.LoopMom[INR],  // KInR
              Var.LoopMom[OUTR], // KOutR
              true);
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
  double Factor = 1.0 / pow(2.0 * PI, D);
  // normalization
  if (LoopNum == 0)
    return 1.0;
  else if (LoopNum == 1) {
    double Tau = Var.Tau[1] - Var.Tau[0];
    double Weight = Fermi.Green(Tau, Var.LoopMom[1], UP, 0, 0);
    Weight *= Fermi.Green(-Tau, Var.LoopMom[1] - Var.LoopMom[0], UP, 0, 0);

    return -SPIN * Weight * Factor;
  }

  // loop order >=2
  dse::polar &P = Polar[LoopNum];
  ver4 &Ver4 = P.Vertex;
  momentum KOutL = Var.LoopMom[1] - Var.LoopMom[0];
  momentum KOutR = Var.LoopMom[2] + Var.LoopMom[0];

  Vertex4(Ver4,
          Var.LoopMom[1], // KInL
          KOutL,          // KOutL
          Var.LoopMom[2], // KInR
          KOutR,          // KOutR
          false);

  // cout << Ver4.Weight[0][DIR] << ", " << Ver4.Weight[0][EX] << endl;
  // cout << Var.LoopMom[0][0] << ", Mom1=" << Var.LoopMom[1].norm() << endl;
  // cout << "DirQ: " << (Var.LoopMom[1] - KOutL).norm() << endl;

  // evaluate all possible G
  EvaluateG(P.G[INL], Var.LoopMom[1]);
  EvaluateG(P.G[OUTL], KOutL);
  EvaluateG(P.G[INR], Var.LoopMom[2]);
  EvaluateG(P.G[OUTR], KOutR);

  int Size = Ver4.T.size();
  P.Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &Gidx = P.Gidx[i];
    double Weight =
        (SPIN * SPIN * Ver4.Weight[i][DIR] + SPIN * Ver4.Weight[i][EX]);

    // attach four G
    for (int j = 0; j < 4; ++j) {
      Weight *= P.G[j][Gidx[j]].Weight;
      // cout << j << "=" << P.G[j][Gidx[j]].Weight << endl;
    }

    P.Weight += Weight;
  }
  // cout << P.Weight << endl;

  return P.Weight * Factor * Factor;
}

double weight::EvaluateSigma(int LoopNum, bool IsFast) {
  // normalization
  double Factor = 1.0 / pow(2.0 * PI, D);
  if (LoopNum == 0) {
    return 1.0;
  } else if (LoopNum == 1) {
    double GWeight =
        Fermi.Green(-EPS, Var.LoopMom[0] + Var.LoopMom[1], UP, 0, 0);
    double VerWeight = VerQTheta.Interaction(Var.LoopMom[1]);
    // cout << GWeight << ", " << VerWeight << endl;
    return GWeight * VerWeight * Factor;
  } else {
    // Sigma with LoopNum>=2
    dse::sigma &Si = Sigma[LoopNum];
    ver4 &Ver4 = Si.Vertex;

    Vertex4(Ver4, Var.LoopMom[0], // KInL
            Var.LoopMom[1],       // KOutL
            Var.LoopMom[1],       // KInR
            Var.LoopMom[0],       // KOutR
            false);

    EvaluateG(Si.G, Var.LoopMom[1]);
    if (IsFast)
      Si.Weight[0] = 0.0;
    else {
      for (auto &w : Si.Weight)
        w = 0.0;
    }
    int Size = Ver4.Weight.size();
    // Size = 1;
    for (int i = 0; i < Size; ++i) {

      auto &Weight = Ver4.Weight[i];
      int Gidx = Si.Gidx[i];
      double w = (Weight[DIR] + Weight[EX] * SPIN) * Si.G[Gidx].Weight;

      if (Si.T[i] != 0)
        w *= 0.5;

      if (IsFast)
        Si.Weight[0] += w * Factor;
      else
        Si.Weight[Si.SigTidx[i]] += w * Factor;
    }
    return Si.Weight[0];
  }
}

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