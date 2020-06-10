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
  double Factor = 1.0 / pow(2.0 * π, D);

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
        ChanWeight[ChanMap[chan]] +=
            W * GWeight * cos(2.0 * π / Para.Beta * dTau);
        // ChanWeight[ChanMap[chan]] += W * GWeight;
      }
    }
  }
}

void vertex4::_EvalUST_CT(const momentum &KInL, const momentum &KOutL,
                          const momentum &KInR, const momentum &KOutR,
                          bool IsFast) {
  if (ChannelCT.size() > 0) {
    double Factor = 1.0 / pow(2.0 * π, D);
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