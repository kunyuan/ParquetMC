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

#define TIND(Shift, LTau, RTau) ((LTau - Shift) * MaxTauNum + RTau - Shift)

void weight::Initialization() {

  // vector<dse::channel> Chan = {dse::T, dse::U, dse::S};
  dse::channel Chan[4] = {dse::I, dse::T, dse::U, dse::S};
  for (int c = 0; c < 4; c++)
    for (int order = 1; order <= Para.Order; order++) {
      vector<dse::channel> chan = {Chan[c]};
      Ver4Root[order][c] =
          VerDiag.Build(Var.LoopMom, order, chan, dse::caltype::PARQUET);
      LOG_INFO(VerDiag.ToString(Ver4Root[order][c]));
    }
}

ver::weightMatrix weight::Evaluate(int LoopNum, int Channel) {
  static ver::weightMatrix Weight;
  if (LoopNum == 0) {
    // normalization
    Weight(DIR) = 1.0;
    Weight(EX) = 0.0;
  } else {
    // if (Channel != dse::T)
    //   return 0.0;

    Weight.SetZero();
    ver4 &Root = Ver4Root[LoopNum][Channel];
    if (Root.Weight.size() != 0) {

      // if (Para.Counter == 12898) {
      //   cout << Root.ID << endl;
      // }

      *Root.LegK[OUTL] = Var.LoopMom[1] - Var.LoopMom[0];
      *Root.LegK[OUTR] = Var.LoopMom[2] + Var.LoopMom[0];

      // if (Channel == dse::S) {
      //   *Root.LegK[INR] = Var.LoopMom[0] - Var.LoopMom[1];
      //   *Root.LegK[OUTR] = Var.LoopMom[0] - Var.LoopMom[2];
      // } else {
      //   *Root.LegK[OUTL] = Var.LoopMom[1] - Var.LoopMom[0];
      //   *Root.LegK[OUTR] = Var.LoopMom[2] + Var.LoopMom[0];
      // }

      Vertex4(Root);

      double Factor = 1.0 / pow(2.0 * PI, D * LoopNum);
      for (auto &w : Root.Weight) {
        Weight(DIR) += w(DIR) * Factor;
        Weight(EX) += w(EX) * Factor;
      }
      // if (LoopNum == 3 && Channel == dse::I) {
      //   cout << "loopnum: " << Root.LoopNum << endl;
      //   cout << "channel: " << Root.Channel[0] << endl;
      //   cout << Weight << endl;
      // }
      // cout << count << endl;
    }
  }
  return Weight;
}

void weight::Ver0(ver4 &Ver4) {
  array<momentum *, 4> &K = Ver4.LegK;
  double &WeightDir = Ver4.Weight[0](DIR); // direct, reducible
  double &WeightEx = Ver4.Weight[0](EX);   // exchange, irreducible
  // Ver4.Weight[0] = 1.0 / Para.Beta;
  if (Ver4.RexpandBare) {
    // bare+quantum correction
    VerQTheta.Interaction(K, 0.0, 1, WeightDir, WeightEx);
  } else {
    // only bare coupling
    VerQTheta.Interaction(K, 0.0, 0, WeightDir, WeightEx);
  }
  return;
}
void weight::Vertex4(dse::ver4 &Ver4) {
  // cout << Ver4.LoopNum << endl;
  if (Ver4.LoopNum == 0) {
    Ver0(Ver4);
  } else {
    for (auto &w : Ver4.Weight)
      w.SetZero();

    ChanUST(Ver4);
    if (Ver4.LoopNum >= 3)
      ChanI(Ver4);
  }
  return;
}

void weight::ChanUST(dse::ver4 &Ver4) {
  double Weight = 0.0, WeightDir = 0.0, WeightEx = 0.0;
  double Ratio;
  array<momentum *, 4> &LegK0 = Ver4.LegK;

  // if (Ver4.ContainProj) {
  // }

  for (auto &bubble : Ver4.Bubble) {
    auto &G = bubble.G;
    const momentum &K0 = *G[0].K;
    int InTL = bubble.InTL;

    for (auto &chan : bubble.Channel)
      if (bubble.IsProjected)
        bubble.ProjFactor[chan] = 0.0;
      else
        bubble.ProjFactor[chan] = 1.0;

    if (bubble.IsProjected && bubble.HasTU) {
      double DirQ = (*LegK0[INL] - *LegK0[OUTL]).norm();
      double ExQ = (*LegK0[INL] - *LegK0[OUTR]).norm();
      if (DirQ < 1.0 * Para.Kf || ExQ < 1.0 * Para.Kf) {
        Ratio = Para.Kf / (*LegK0[INL]).norm();
        *bubble.LegK[T][INL] = *LegK0[INL] * Ratio;
        Ratio = Para.Kf / (*LegK0[INR]).norm();
        *bubble.LegK[T][INR] = *LegK0[INR] * Ratio;
        if (DirQ < 1.0 * Para.Kf) {
          // *bubble.LegK[T][OUTL] = *bubble.LegK[T][INL];
          // *bubble.LegK[T][OUTR] = *bubble.LegK[T][INR];
          // double x=
          bubble.ProjFactor[T] = exp(-DirQ * DirQ / 0.1);
          // if (DirQ < EPS)
          //   bubble.ProjFactor[T] = 1.0;
        }
        if (ExQ < 1.0 * Para.Kf) {
          // *bubble.LegK[U][OUTL] = *bubble.LegK[T][INR];
          // *bubble.LegK[U][OUTR] = *bubble.LegK[T][INL];
          // if (ExQ < EPS)
          // bubble.ProjFactor[U] = 1.0;
          bubble.ProjFactor[U] = exp(-ExQ * ExQ / 0.1);
        }
      }
    }

    if (bubble.IsProjected && bubble.HasS) {
      // double InL = (*LegK0[INL]).norm();
      // double OutL = (*LegK0[OUTL]).norm();
      // double InR = (*LegK0[INR]).norm();
      // double OutR = (*LegK0[OUTR]).norm();

      // double DirQ = (*LegK0[INL] - *LegK0[OUTL]).norm();
      // if (DirQ < 1.0 * Para.Kf) {
      //   Ratio = Para.Kf / (*LegK0[INL]).norm();
      //   *bubble.LegK[S][INL] = *LegK0[INL] * Ratio;
      //   Ratio = Para.Kf / (*LegK0[INR]).norm();
      //   *bubble.LegK[S][INR] = *LegK0[INR] * Ratio;
      //   *bubble.LegK[S][OUTL] = *bubble.LegK[S][INL];
      //   *bubble.LegK[S][OUTR] = *bubble.LegK[S][INR];
      //   bubble.ProjFactor[S] = exp(-DirQ * DirQ / 0.1);
      // }
      // if ((InL < 1.1 * Para.Kf && InL > 0.9 * Para.Kf) &&
      //     (OutL < 1.1 * Para.Kf && OutL > 0.9 * Para.Kf) &&
      //     (InR < 1.1 * Para.Kf && InR > 0.9 * Para.Kf) &&
      //     (OutR < 1.1 * Para.Kf && OutR > 0.9 * Para.Kf)) {
      //   Ratio = Para.Kf / (*LegK0[INL]).norm();
      //   *bubble.LegK[S][INL] = *LegK0[INL] * Ratio;
      //   Ratio = Para.Kf / (*LegK0[INR]).norm();
      //   *bubble.LegK[S][INR] = *LegK0[INR] * Ratio;
      //   Ratio = Para.Kf / (*LegK0[OUTL]).norm();
      //   *bubble.LegK[S][OUTL] = *LegK0[OUTL] * Ratio;
      //   Ratio = Para.Kf / (*LegK0[OUTR]).norm();
      //   *bubble.LegK[S][OUTR] = *LegK0[OUTR] * Ratio;
      //   bubble.ProjFactor[S] = 1.0;
      // }
      //   double InQ = (*LegK0[INL] + *LegK0[INR]).norm();
      //   if (InQ < 1.0 * Para.Kf) {
      //     Ratio = Para.Kf / (*LegK0[INL]).norm();
      //     *bubble.LegK[S][INL] = *LegK0[INL] * Ratio;
      //     Ratio = Para.Kf / (*LegK0[OUTL]).norm();
      //     *bubble.LegK[S][OUTL] = *LegK0[OUTL] * Ratio;
      //     *bubble.LegK[S][INR] = *bubble.LegK[S][INL] * (-1.0);
      //     *bubble.LegK[S][OUTR] = *bubble.LegK[S][OUTL] * (-1.0);
      //     bubble.ProjFactor[S] = exp(-InQ * InQ / 0.1);
      //   }
    }

    for (auto &chan : bubble.Channel) {
      array<momentum *, 4> &LegK = bubble.LegK[chan];
      if (chan == T)
        *G[T].K = *LegK[OUTL] + K0 - *LegK[INL];
      else if (chan == U)
        *G[U].K = *LegK[OUTR] + K0 - *LegK[INL];
      else if (chan == S)
        *G[S].K = *LegK[INL] + *LegK[INR] - K0;
    }

    for (int lt = InTL; lt < InTL + Ver4.TauNum - 1; ++lt)
      for (int rt = InTL + 1; rt < InTL + Ver4.TauNum; ++rt) {
        double dTau = Var.Tau[rt] - Var.Tau[lt];
        G[0](lt, rt) = Fermi.Green(dTau, K0, UP, 0, Var.CurrScale);
        for (auto &chan : bubble.Channel) {
          // if (chan > 3)
          //   ABORT("too many chan " << chan);
          if (abs(bubble.ProjFactor[chan]) > EPS)
            if (chan == S)
              // LVer to RVer
              G[S](lt, rt) = Fermi.Green(dTau, *G[S].K, UP, 0, Var.CurrScale);
            else
              // RVer to LVer
              G[chan](rt, lt) =
                  Fermi.Green(-dTau, *G[chan].K, UP, 0, Var.CurrScale);
        }
      }

    // for vertex4 with one or more loops
    for (auto &pair : bubble.Pair) {
      if (abs(bubble.ProjFactor[pair.Channel]) < EPS)
        continue;
      ver4 &LVer = pair.LVer;
      ver4 &RVer = pair.RVer;
      Vertex4(LVer);
      Vertex4(RVer);

      for (auto &map : pair.Map) {
        Weight = pair.SymFactor * bubble.ProjFactor[pair.Channel];
        Weight *= G[0](map.G0T) * G[pair.Channel](map.GT);
        auto &CWeight = Ver4.Weight[map.Tidx];
        auto &LWeight = LVer.Weight[map.LVerTidx];
        auto &RWeight = RVer.Weight[map.RVerTidx];
        if (pair.Channel == T) {
          CWeight(DIR) += Weight * (LWeight(DIR) * RWeight(DIR) * SpinIndex +
                                    LWeight(DIR) * RWeight(EX) +
                                    LWeight(EX) * RWeight(DIR));
          CWeight(EX) += Weight * LWeight(EX) * RWeight(EX);
        } else if (pair.Channel == U) {
          CWeight(EX) += Weight * (LWeight(DIR) * RWeight(DIR) * SpinIndex +
                                   LWeight(DIR) * RWeight(EX) +
                                   LWeight(EX) * RWeight(DIR));
          CWeight(DIR) += Weight * LWeight(EX) * RWeight(EX);
        } else if (pair.Channel == S) {
          CWeight(DIR) += Weight * (LWeight(DIR) * RWeight(DIR) +
                                    LWeight(EX) * RWeight(EX));

          CWeight(EX) += Weight * (LWeight(DIR) * RWeight(EX) +
                                   LWeight(EX) * RWeight(DIR));
        }
      }
    }
  }
}

void weight::ChanI(dse::ver4 &Ver4) {
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
  //     g.Weight = Fermi.Green(Var.Tau[g.OutT] - Var.Tau[g.InT], *(g.K), UP, 0,
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
