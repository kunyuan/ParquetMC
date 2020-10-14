#include "vertex4.h"
#include <algorithm>
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
  _EvalI(KInL, KOutL, KInR, KOutR, IsFast);
  return;
}

void vertex4::_EvalBare(const momentum &KInL, const momentum &KOutL,
                        const momentum &KInR, const momentum &KOutR) {

  // only bare coupling
  if (DiagType == POLAR)
    Weight[0] =
        Prop.Interaction(KInL, KOutL, KInR, KOutR, Var.LoopMom[0].norm());
  else if ((DiagType == GAMMA) && IsProper)
    Weight[0] = Prop.Interaction(KInL, KOutL, KInR, KOutR,
                                 (Var.LoopMom[INL] - Var.LoopMom[OUTL]).norm());
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
        if (IsF) {
          double dTau = Var.Tau[Tpair[t][INL]] + Var.Tau[Tpair[t][OUTL]];
          dTau += -Var.Tau[Tpair[t][INR]] - Var.Tau[Tpair[t][OUTR]];
          ChanWeight[ChanMap[chan]] += W * GWeight * cos(PI / Para.Beta * dTau);
        } else {
          double dTau = Var.Tau[Tpair[t][INL]] - Var.Tau[Tpair[t][OUTL]];
          dTau += Var.Tau[Tpair[t][INR]] - Var.Tau[Tpair[t][OUTR]];
          ChanWeight[ChanMap[chan]] += W * GWeight * cos(PI / Para.Beta * dTau);
        }
      }
    }
  }
}

// void vertex4::_EvalUST_CT(const momentum &KInL, const momentum &KOutL,
//                           const momentum &KInR, const momentum &KOutR,
//                           bool IsFast) {
//   if (ChannelCT.size() > 0) {
//     double Factor = 1.0 / pow(2.0 * PI, D);
//     // cout << weight << endl;
//     // ProjFactor = SymFactor[chan] * Factor;
//     // cout << Tpair[0][0] << ", " << Tpair[0][1] << Tpair[0][2] <<
//     Tpair[0][3]
//     //      << endl;
//     for (auto &c : ChannelCT) {
//       double weight;
//       if (c == TC) {
//         if (DiagType == POLAR)
//           weight = pow(Prop.Interaction(KInL - KOutL, 0,
//           Var.LoopMom[0].norm()),
//                        LoopNum() + 1);
//         else if (IsProper && DiagType == GAMMA)
//           weight = pow(
//               Prop.Interaction(KInL - KOutL, 0,
//                                (Var.LoopMom[INL] -
//                                Var.LoopMom[OUTL]).norm()),
//               LoopNum() + 1);
//         else
//           weight = pow(Prop.Interaction(KInL - KOutL, 0), LoopNum() + 1);

//         for (int o = LoopIdx; o < LoopIdx + LoopNum(); o++) {
//           weight *= Prop.CounterBubble(Var.LoopMom[o]) * Factor * SPIN;
//         }
//         if (IsFast && Level == 0)
//           ChanWeight[T][DIR] += weight * SymFactor[TC];
//         else
//           // counter-term Tpair and the weight are always the first element
//           Weight[0][DIR] += weight * SymFactor[TC];
//       } else if (c == UC) {
//         if (DiagType == POLAR)
//           weight = pow(Prop.Interaction(KInL - KOutR, 0,
//           Var.LoopMom[0].norm()),
//                        LoopNum() + 1);
//         else if (IsProper && DiagType == GAMMA)
//           weight = pow(
//               Prop.Interaction(KInL - KOutR, 0,
//                                (Var.LoopMom[INL] -
//                                Var.LoopMom[OUTL]).norm()),
//               LoopNum() + 1);
//         else
//           weight = pow(Prop.Interaction(KInL - KOutR, 0), LoopNum() + 1);

//         // double weight =
//         //     pow(Prop.Interaction(KInL - KOutR, 0, Var.LoopMom[0].norm()),
//         //         LoopNum() + 1);

//         for (int o = LoopIdx; o < LoopIdx + LoopNum(); o++) {
//           weight *= Prop.CounterBubble(Var.LoopMom[o]) * Factor * SPIN;
//         }
//         if (IsFast && Level == 0)
//           ChanWeight[U][EX] += weight * SymFactor[UC];
//         else
//           // counter-term Tpair and the weight are always the first element
//           Weight[0][EX] += weight * SymFactor[UC];
//       }
//     }
//   }
// }

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
      double wd, we, bubbles;
      double K = (KInL - KOutL).norm();
      if (c == TC) {
        if (DiagType == POLAR) {
          bool isproper = IsEqual(K, Var.LoopMom[0].norm());
          // if IsProper=true, then only one-interaction irreducible diagrams
          // are allowed;
          we = SPIN * pow(Prop.Rm(0.0, K, true, isproper), LoopNum() + 1);
          wd = pow(Prop.Rp(0.0, K, true, isproper), LoopNum() + 1) - we / SPIN;
        }
        //  else if (IsProper && DiagType == GAMMA)
        //   weight = pow(
        //       Prop.Interaction(KInL - KOutL, 0,
        //                        (Var.LoopMom[INL] -
        //                        Var.LoopMom[OUTL]).norm()),
        //       LoopNum() + 1);
        // else
        //   weight = pow(Prop.Interaction(KInL - KOutL, 0), LoopNum() + 1);

        bubbles = 1.0;
        for (int o = LoopIdx; o < LoopIdx + LoopNum(); o++) {
          bubbles *= Prop.CounterBubble(Var.LoopMom[o]) * Factor * SPIN;
        }
        // bubbles *= pow(-1.0, LoopNum() + 1);
        if (IsFast && Level == 0) {
          ChanWeight[T][DIR] += wd * bubbles * SymFactor[TC];
          ChanWeight[T][EX] += we * bubbles * SymFactor[TC];
        } else {
          // counter-term Tpair and the weight are always the first element
          Weight[0][DIR] += wd * bubbles * SymFactor[TC];
          Weight[0][EX] += we * bubbles * SymFactor[TC];
        }
      } else if (c == UC) {
        double weight;
        if (DiagType == POLAR)
          weight = pow(Prop.Interaction(KInL - KOutR, 0, Var.LoopMom[0].norm()),
                       LoopNum() + 1);
        else if (IsProper && DiagType == GAMMA)
          weight = pow(
              Prop.Interaction(KInL - KOutR, 0,
                               (Var.LoopMom[INL] - Var.LoopMom[OUTL]).norm()),
              LoopNum() + 1);
        else
          weight = pow(Prop.Interaction(KInL - KOutR, 0), LoopNum() + 1);

        // double weight =
        //     pow(Prop.Interaction(KInL - KOutR, 0, Var.LoopMom[0].norm()),
        //         LoopNum() + 1);

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
                     const momentum &KInR, const momentum &KOutR, bool IsFast) {
  return;
  // only loop 3 has envolope diagrams
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

  double G0 = Prop.Green(T2 - T0, K0, UP);
  double G1 = Prop.Green(T1 - T0, K1, UP);
  double G2 = Prop.Green(T1 - T2, K2, UP);
  double G3 = Prop.Green(T0 - T3, K3, UP);
  double G4 = Prop.Green(T3 - T1, K4, UP);
  double G5 = Prop.Green(T3 - T2, K5, UP);
  double G6 = Prop.Green(T3 - T1, K6, UP);
  double G7 = Prop.Green(T2 - T3, K7, UP);
  double G8 = Prop.Green(T2 - T3, K8, UP);

  verWeight v0 = Prop.Interaction(KInL, K1, K3, K0);
  verWeight v1 = Prop.Interaction(K1, KOutL, K2, K4);
  verWeight v2 = Prop.Interaction(K1, KOutR, K2, K6);
  verWeight v3 = Prop.Interaction(K0, K2, KInR, K5);
  verWeight v4 = Prop.Interaction(K0, K2, K7, KOutR);
  verWeight v5 = Prop.Interaction(K0, K2, K8, KOutL);
  verWeight v6 = Prop.Interaction(K4, K3, K5, KOutR);
  verWeight v7 = Prop.Interaction(K6, K3, K5, KOutL);
  verWeight v8 = Prop.Interaction(K4, K3, KInR, K7);
  verWeight v9 = Prop.Interaction(K6, K3, KInR, K8);

  double Factor = 1.0 / pow(2.0 * PI, 3.0 * D);
  verWeight VerWeight;

  // Diagram 1
  double GWeight = G0 * G1 * G2 * G3 * G4 * G5;
  VerWeight = (-1.0) * Factor * GWeight * _Envolope12(v0, v1, v3, v6, false);
  if (IsFast && Level == 0)
    ChanWeight[I] += VerWeight * cos(PI / Para.Beta * (T0 - T1 + T2 - T3));
  else
    Weight[_EnvolpeVerIdx[0]] += VerWeight;

  // Diagram 2
  GWeight = G0 * G1 * G2 * G3 * G5 * G6;
  VerWeight = (+1.0) * Factor * GWeight * _Envolope12(v0, v2, v3, v7, true);
  if (IsFast && Level == 0)
    ChanWeight[I] += VerWeight * cos(PI / Para.Beta * (T0 - T3 + T2 - T1));
  else
    Weight[_EnvolpeVerIdx[1]] += VerWeight;

  // Diagram 3
  GWeight = G0 * G1 * G2 * G3 * G4 * G7;
  VerWeight = (-1.0) * Factor * GWeight * _Envolope34(v0, v1, v4, v8, false);
  if (IsFast && Level == 0)
    ChanWeight[I] += VerWeight * cos(PI / Para.Beta * (T0 - T1 + T3 - T2));
  else
    Weight[_EnvolpeVerIdx[2]] += VerWeight;

  // Diagram 4
  GWeight = G0 * G1 * G2 * G3 * G6 * G8;
  VerWeight = (+1.0) * Factor * GWeight * _Envolope34(v0, v2, v5, v9, true);
  if (IsFast && Level == 0)
    ChanWeight[I] += VerWeight * cos(PI / Para.Beta * (T0 - T2 + T3 - T1));
  else
    Weight[_EnvolpeVerIdx[3]] += VerWeight;
}

verWeight vertex4::_Envolope12(const verWeight &iL, const verWeight &oL,
                               const verWeight &iR, const verWeight &oR,
                               bool isexchange) {
  verWeight W = {0.0, 0.0};
  int S = SPIN;
  double iLd = iL[DIR], iLe = iL[EX];
  double oLd = oL[DIR], oLe = oL[EX];
  double iRd = iR[DIR], iRe = iR[EX];
  double oRd = oR[DIR], oRe = oR[EX];
  // diagram 11, 12, 13, 14
  W[DIR] += iLd * oLd * (S * (iRd * oRd + iRe * oRe) + iRe * oRd + iRd * oRe);
  // diagram 21,31, 41
  W[DIR] += (iLd * oLe + iLe * oLd + S * iLe * oLe) * iRd * oRd;
  // diagram 32, 43
  W[DIR] += iLe * oLd * iRe * oRe + iLe * oLe * iRd * oRe;

  // diagram 22, 23, 24
  W[EX] += iLd * oLe * (S * iRe * oRe + iRd * oRe + iRe * oRd);
  // diagram 33, 34
  W[EX] += iLe * oLd * (iRd * oRe + iRe * oRd);
  // diagram 42, 44
  W[EX] += iLe * oLe * (iRe * oRe + S * iRe * oRd);
  if (isexchange) {
    swap(W[DIR], W[EX]);
    // double temp = W[DIR];
    // W[DIR] = W[EX];
    // W[EX] = temp;
  }
  // swap(W[DIR], W[EX]);
  return W;
}

verWeight vertex4::_Envolope34(const verWeight &iL, const verWeight &oL,
                               const verWeight &iR, const verWeight &oR,
                               bool isexchange) {
  verWeight W = {0.0, 0.0};
  int S = SPIN;
  double iLd = iL[DIR], iLe = iL[EX];
  double oLd = oL[DIR], oLe = oL[EX];
  double iRd = iR[DIR], iRe = iR[EX];
  double oRd = oR[DIR], oRe = oR[EX];

  // diagram 11, 12, 13, 14
  W[DIR] += iLd * oLd * (S * (iRd * oRd + iRe * oRe) + iRe * oRd + iRd * oRe);
  // diagram 21,31, 41
  W[DIR] += (iLd * oLe + iLe * oLd + S * iLe * oLe) * iRd * oRd;
  // diagram 22, 43
  W[DIR] += iLd * oLe * iRe * oRe + iLe * oLe * iRd * oRe;

  // diagram 23, 24
  W[EX] += iLd * oLe * (iRd * oRe + iRe * oRd);
  // diagram 32, 33, 34
  W[EX] += iLe * oLd * (S * iRe * oRe + iRd * oRe + iRe * oRd);
  // diagram 42, 44
  W[EX] += iLe * oLe * (iRe * oRe + S * iRe * oRd);
  if (isexchange) {
    swap(W[DIR], W[EX]);
    // double temp = W[DIR];
    // W[DIR] = W[EX];
    // W[EX] = temp;
  }
  return W;
}