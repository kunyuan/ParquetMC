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
    Ver4Root[order] =
        VerDiag.Vertex(0, order, 4, 0, Chan, ExtLegK, RIGHT, false);
    LOG_INFO(VerDiag.ToString(Ver4Root[order]));
  }
}

ver::weightMatrix weight::Evaluate(int LoopNum, int Channel) {
  ver::weightMatrix Weight;
  if (LoopNum == 0) {
    // normalization
    Weight[EX] = 0.0;
    Weight[DIR] = 1.0;
  } else {
    // if (Channel != dse::T)
    //   return 0.0;
    Weight.SetZero();

    ver4 &Root = Ver4Root[LoopNum];
    if (Root.Weight.size() != 0) {

      for (int c = 0; c < 4; ++c)
        Root.ChanWeight[c].SetZero();

      Vertex4(Root);

      if (Channel <= 4) {
        double Factor = 1.0 / pow(2.0 * PI, D * LoopNum);

        //////// Measure Scattering amplitude ////////////////////
        for (auto &w : Root.Weight) {
          Weight[DIR] += w[DIR] * Factor;
          Weight[EX] += w[EX] * Factor;
        }
      }
    }
  }
  return Weight;
}

void weight::Ver0(ver4 &Ver4) {
  const auto &LegK = Ver4.LegK;
  double &WeightDir = Ver4.Weight[0][DIR]; // direct, reducible
  double &WeightEx = Ver4.Weight[0][EX];   // exchange, irreducible
  // Ver4.Weight[0] = 1.0 / Para.Beta;
  if (Para.Type == BARE || (Para.Type == LEFT && Para.Type == PARQUET)) {
    // only bare coupling
    VerQTheta.Interaction(LegK, 0.0, false, Ver4.InBox, WeightDir, WeightEx);
  } else {
    // bare+quantum correction
    VerQTheta.Interaction(LegK, 0.0, true, Ver4.InBox, WeightDir, WeightEx);
  }
  return;
}
void weight::Vertex4(ver4 &Ver4) {
  // cout << Ver4.LoopNum << endl;
  if (Ver4.LoopNum == 0) {
    Ver0(Ver4);
  } else {
    for (auto &w : Ver4.Weight)
      w.SetZero();

    // calculate four channels for the root vertex
    if (Ver4.Level == 0)
      for (auto &w : Ver4.ChanWeight)
        w.SetZero();

    ChanUST(Ver4);
    if (Ver4.LoopNum >= 3)
      ChanI(Ver4);
  }
  return;
}

void weight::ChanUST(ver4 &Ver4) {
  double Weight = 0.0;
  double Ratio, dTau;
  double DirW, ExW;
  const auto &LegK = Ver4.LegK;
  auto &G = Ver4.G;

  for (auto chan : Ver4.Channel) {
    if (chan < 4) {

      if (chan == 0)
        Ver4.K[0] = Var.LoopMom[Ver4.Loopidx];
      else if (chan == T)
        Ver4.K[T] = *LegK[OUTL] + Ver4.K[0] - *LegK[INL];
      else if (chan == U)
        Ver4.K[U] = *LegK[OUTR] + Ver4.K[0] - *LegK[INL];
      else if (chan == S)
        Ver4.K[S] = *LegK[INL] + *LegK[INR] - Ver4.K[0];

      for (auto &g : G[chan]) {
        dTau = Var.Tau[g.T[OUT]] - Var.Tau[g.T[IN]];
        g.Weight = Fermi.Green(dTau, Ver4.K[chan], UP, 0, Var.CurrScale);
      }
    }
  }
  ///////////// Check if the projected counter-terms exist or not ///////

  // for vertex4 with one or more loops
  for (auto &pair : Ver4.Pair) {
    ver4 &LVer = pair.LVer;
    ver4 &RVer = pair.RVer;
    Vertex4(LVer);
    Vertex4(RVer);

    for (auto &map : pair.Map) {
      Weight = SymFactor[pair.Channel];
      Weight *= G[0][map.G0idx].Weight * G[pair.Channel][map.Gidx].Weight;

      auto &Lw = LVer.Weight[map.LVerTidx];
      auto &Rw = RVer.Weight[map.RVerTidx];

      switch (pair.Channel) {
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
      default:
        ABORT("Channel does not exist!");
        break;
      }

      Ver4.Weight[map.Tidx][DIR] += DirW * Weight;
      Ver4.Weight[map.Tidx][EX] += ExW * Weight;

      if (Ver4.Level == 0) {
        Ver4.ChanWeight[pair.Channel][DIR] += DirW * Weight;
        Ver4.ChanWeight[pair.Channel][EX] += ExW * Weight;
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
