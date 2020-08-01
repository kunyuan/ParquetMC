#include "vertex4.h"
#include <iostream>

using namespace std;
using namespace diag;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void vertex4::Evaluate(const momentum &KInL, const momentum &KOutL,
                       const momentum &KInR, const momentum &KOutR,
                       bool IsFast) {

  if (Order == 0) {
    _EvalBare(KInL, KOutL, KInR, KOutR);
    return;
  }
  if (Level > 0 || IsFast == false) {
    for (auto &w : Weight)
      w.setZero();
  } else {
    for (auto &w : ChanWeight)
      w.setZero();
  }

  _EvalUST(KInL, KOutL, KInR, KOutR, IsFast);
  _EvalUST_CT(KInL, KOutL, KInR, KOutR, IsFast);
  // if (Ver4.LoopNum >= 3)
  //   _EvalI(KInL, KOutL, KInR, KOutR, IsFast);
  return;
}

void vertex4::_EvalBare(const momentum &KInL, const momentum &KOutL,
                        const momentum &KInR, const momentum &KOutR) {

  // only bare coupling
  if (DiagType == POLAR)
    Weight[0] =
        Prop.Interaction(KInL, KOutL, KInR, KOutR, Var.LoopMom[0].norm());
  else
    Weight[0] = Prop.Interaction(KInL, KOutL, KInR, KOutR);
  return;
}

void vertex4::_EvalUST(const momentum &KInL, const momentum &KOutL,
                       const momentum &KInR, const momentum &KOutR,
                       bool IsFast) {
  // calculate K table
  G[0].K = Var.LoopMom[LoopIdx];
  G[0].Evaluate();

  for (auto chan : Channel) {
    if (chan == T) {
      G[T].K = KOutL + G[0].K - KInL;
      G[T].Evaluate();
    } else if (chan == U) {
      G[U].K = KOutR + G[0].K - KInL;
      G[U].Evaluate();
    } else if (chan == S) {
      G[S].K = KInL + KInR - G[0].K;
      G[S].Evaluate();
    }
  }

  // for vertex4 with one or more loops
  verWeight W;
  double GWeight, ProjFactor;
  double Factor = 1.0 / pow(2.0 * PI, D);

  for (auto &b : _UST) {
    // cout << "before, " << b.Channel << ", " << b.Map.size() << endl;
    channel chan = b.Channel;
    ProjFactor = SymFactor[chan] * Factor;

    // cout << _UST.size() << endl;
    // cout << "before: " << b.Channel << ", " << b.Map.size() << endl;
    if (chan == T) {
      b.LVer.Evaluate(KInL, KOutL, G[T].K, G[0].K);
      b.RVer.Evaluate(G[0].K, G[T].K, KInR, KOutR);
    } else if (chan == U) {
      b.LVer.Evaluate(KInL, KOutR, G[U].K, G[0].K);
      b.RVer.Evaluate(G[0].K, G[U].K, KInR, KOutL);
    } else if (chan == S) {
      b.LVer.Evaluate(KInL, G[S].K, KInR, G[0].K);
      b.RVer.Evaluate(G[0].K, KOutL, G[S].K, KOutR);
    }
    // cout << "middle: " << b.Channel << ", " << b.Map.size() << endl;

    auto &LVerW = b.LVer.Weight;
    auto &RVerW = b.RVer.Weight;
    // cout << "after1: " << b.Channel << endl;
    // cout << "after2: " << b.Channel << ", " << b.Map.size() << endl << endl;

    for (auto &map : b.Map) {
      GWeight = ProjFactor * G[0][map[G0T]] * G[chan][map[GXT]];

      auto &Lw = LVerW[map[LVERT]];
      auto &Rw = RVerW[map[RVERT]];

      // cout << "g " << G[0][map[G0T]] << ", " << G[chan][map[GXT]] << endl;
      // cout << "Order: " << Order << ", c=" << chan << ", " << GWeight << ",
      // ["
      //      << Lw[DIR] << ", " << Lw[EX] << "], [" << Rw[DIR] << ", " <<
      //      Rw[EX]
      //      << "]" << endl;

      if (chan == T) {
        W[DIR] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        W[EX] = Lw[EX] * Rw[EX];
      } else if (chan == U) {
        W[EX] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        W[DIR] = Lw[EX] * Rw[EX];
      } else if (chan == S) {
        // see the note "code convention"
        W[DIR] = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
        W[EX] = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX];
      }

      if (Level > 0 || IsFast == false)
        Weight[map[VERT]] += W * GWeight;
      else {
        int t = map[VERT];
        double dTau = Var.Tau[Tpair[t][INL]] - Var.Tau[Tpair[t][OUTL]];
        dTau += Var.Tau[Tpair[t][INR]] - Var.Tau[Tpair[t][OUTR]];
        ChanWeight[ChanMap[chan]] += W * GWeight * cos(PI / Para.Beta * dTau);
        // ChanWeight[ChanMap[chan]] += W * GWeight;
      }
    }
  }
}

void vertex4::_EvalUST_CT(const momentum &KInL, const momentum &KOutL,
                          const momentum &KInR, const momentum &KOutR,
                          bool IsFast) {
  if (ChannelCT.size() > 0) {
    double Factor = 1.0 / pow(2.0 * PI, D);
    // cout << weight << endl;
    // ProjFactor = SymFactor[chan] * Factor;
    // cout << Tpair[0][0] << ", " << Tpair[0][1] << Tpair[0][2] << Tpair[0][3]
    //      << endl;
    for (auto &c : ChannelCT) {
      if (c == TC) {
        double weight =
            pow(Prop.Interaction(KInL - KOutL, 0, Var.LoopMom[0].norm()),
                LoopNum() + 1);
        for (int o = LoopIdx; o < LoopIdx + LoopNum(); o++) {
          weight *= Prop.CounterBubble(Var.LoopMom[o]) * Factor * SPIN;
        }
        if (IsFast && Level == 0)
          ChanWeight[T][DIR] += weight * SymFactor[TC];
        else
          // counter-term Tpair and the weight are always the first element
          Weight[0][DIR] += weight * SymFactor[TC];
      } else if (c == UC) {
        double weight =
            pow(Prop.Interaction(KInL - KOutR, 0, Var.LoopMom[0].norm()),
                LoopNum() + 1);
        for (int o = LoopIdx; o < LoopIdx + LoopNum(); o++) {
          weight *= Prop.CounterBubble(Var.LoopMom[o]) * Factor * SPIN;
        }
        if (IsFast && Level == 0)
          ChanWeight[U][EX] += weight * SymFactor[UC];
        else
          // counter-term Tpair and the weight are always the first element
          Weight[0][EX] += weight * SymFactor[UC];
      }
    }
  }
}

void vertex4::_EvalI(const momentum &KInL, const momentum &KOutL,
                     const momentum &KInR, const momentum &KOutR,
                     bool IsFast = false) {
  if (Order != 3)
    return;

  momentum &K0 = Var.LoopMom[LoopIdx];
  momentum &K1 = Var.LoopMom[LoopIdx + 1];
  momentum &K2 = Var.LoopMom[LoopIdx + 2];
  momentum K3 = K0 + K1 - KInL;
  momentum K4 = K1 + K2 - KOutL;
  momentum K5 = K0 + KInR - K2;
  momentum K6 = K1 + K2 - KOutR;
  momentum K7 = K2 + KOutR - K0;
  momentum K8 = K2 + KOutL - K0;

  double T0 = Var.Tau[Tidx];
  double T1 = Var.Tau[Tidx + 1];
  double T2 = Var.Tau[Tidx + 2];
  double T3 = Var.Tau[Tidx + 3];

  double G1 = Prop.Green(T1 - T0, K1, UP);
  double G2 = Prop.Green(T1 - T2, K2, UP);
  double G3 = Prop.Green(T0 - T3, K3, UP);
  double G4 = Prop.Green(T3 - T1, K4, UP);
  double G5 = Prop.Green(T3 - T2, K5, UP);
  double G6 = Prop.Green(T3 - T1, K6, UP);
  double G7 = Prop.Green(T2 - T3, K7, UP);
  double G8 = Prop.Green(T2 - T3, K8, UP);

  verWeight ver0 = Prop.Interaction(KInL, K1, K3, K0);
  verWeight ver1 = Prop.Interaction(K1, KOutL, K2, K4);
  verWeight ver2 = Prop.Interaction(K1, KOutR, K2, K6);
  verWeight ver3 = Prop.Interaction(K0, K2, KInR, K5);
  verWeight ver4 = Prop.Interaction(K0, K2, K7, KOutR);
  verWeight ver5 = Prop.Interaction(K0, K2, K8, KOutL);
  verWeight ver6 = Prop.Interaction(K4, K3, K5, KOutR);
  verWeight ver7 = Prop.Interaction(K6, K3, K5, KOutL);
  verWeight ver8 = Prop.Interaction(K4, K3, KInR, K7);
  verWeight ver9 = Prop.Interaction(K6, K3, KInR, K8);

  double Weight = 0.0;
  double ComWeight = 0.0;
  for (auto &map : Env.Map) {
    auto &SubVer = Env.Ver;
    auto &GT = map.GT;
    auto &G = Env.G;
    ComWeight = G[0].Weight * G[1].Weight * G[2].Weight * G[3].Weight;
    // cout << "G: " << ComWeight << endl;
    ComWeight *= SubVer[0].Weight[map.LDVerTidx];
    // cout << "Ver: " << SubVer[0].Weight[map.LDVerT] << endl;
    // cout << "T: " << map.LDVerT << endl;

    Weight = Env.SymFactor[0] * ComWeight;
    Weight *= SubVer[1].Weight[map.LUVerTidx];
    Weight *= SubVer[3].Weight[map.RDVerTidx];
    Weight *= SubVer[6].Weight[map.RUVerTidx];
    Weight *= G[4].Weight * G[5].Weight;
    Ver4.Weight[map.Tidx[0]] += Weight;

    Weight = Env.SymFactor[1] * ComWeight;
    Weight *= SubVer[2].Weight[map.LUVerTidx];
    Weight *= SubVer[3].Weight[map.RDVerTidx];
    Weight *= SubVer[7].Weight[map.RUVerTidx];
    Weight *= G[6].Weight * G[5].Weight;
    Ver4.Weight[map.Tidx[1]] += Weight;
    // cout << Weight << endl;

    Weight = Env.SymFactor[2] * ComWeight;
    Weight *= SubVer[1].Weight[map.LUVerTidx];
    Weight *= SubVer[4].Weight[map.RDVerTidx];
    Weight *= SubVer[8].Weight[map.RUVerTidx];
    Weight *= G[4].Weight * G[7].Weight;
    Ver4.Weight[map.Tidx[2]] += Weight;
    // cout << Weight << endl;

    Weight = Env.SymFactor[3] * ComWeight;
    Weight *= SubVer[2].Weight[map.LUVerTidx];
    Weight *= SubVer[5].Weight[map.RDVerTidx];
    Weight *= SubVer[9].Weight[map.RUVerTidx];
    Weight *= G[6].Weight * G[8].Weight;
    Ver4.Weight[map.Tidx[3]] += Weight;
    // cout << Weight << endl;

    // if (map.LDVerT == 0 && map.LUVerT == 0 && map.RDVerT == 0 &&
    //     map.RUVerT == 0) {
    // cout << "Com: " << ComWeight << endl;
    // cout << "G[4]: " << G[4](GT[4]) << endl;
    // cout << "G[5]: " << G[5](GT[5]) << endl;
    // cout << SubVer[1].Weight[map.LUVerT] << endl;
    // cout << SubVer[3].Weight[map.RDVerT] << endl;
    // cout << SubVer[6].Weight[map.RUVerT] << endl;
    // cout << "First: " << Weight << endl;
    // }
  }
}