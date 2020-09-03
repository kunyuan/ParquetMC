#define FMT_HEADER_ONLY
#include "diagram.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void sigma::Build(int order) {
  ASSERT_ALLWAYS(order >= 0, "Sigma LoopNum must be larger than 0!");
  Order = order;
  if (order <= 1)
    return;
  ASSERT(TauNum() < MaxTauNum - 1, "MaxTauNum is too small!");

  vector<channel> FULL = {I, T, U, S, TC, UC};
  vector<channel> F = {TC, UC};
  // the bare interaction is automatically included

  for (int ol = 0; ol < LoopNum() - 1; ol++) {
    verPair bub;

    ////////////////////   Right SubVer  ///////////////////
    int lvl = 0;
    int oR = LoopNum() - 2 - ol;
    int LInTL = 0;
    int RInTL = MaxTauNum - 1;
    int Llopidx = 3; // ExtK: 0, G1: 1, G2: 2
    int Rlopidx = 3 + ol;

    bub.LVer.Build(lvl, ol, Llopidx, LInTL, FULL, LEFT);
    bub.RVer.Build(lvl, oR, Rlopidx, RInTL, F, RIGHT);

    for (int ol = 0; ol < bub.LVer.Tpair.size(); ++ol) {
      //   cout << ol << endl;
      // Assume the RVer is equal-time!
      auto &t = bub.LVer.Tpair[ol];

      int G1idx = G1.AddTidxPair({t[OUTL], RInTL});
      int G2idx = G2.AddTidxPair({t[OUTR], RInTL});
      int G3idx = G3.AddTidxPair({RInTL, t[INR]});

      bub.Map.push_back({ol, G1idx, G2idx, G3idx});
    }

    // cout << Map.size() << endl;
    Bubble.push_back(bub);
  }
}

double sigma::Evaluate() {
  // normalization
  double Factor = 1.0 / pow(2.0 * PI, D);
  if (Order == 0) {
    return 1.0;
  } else if (Order == 1) {
    double GWeight = Prop.Green(-EPS, Var.LoopMom[0] + Var.LoopMom[1], UP, 0);
    double VerWeight = Prop.Interaction(Var.LoopMom[1], 0);
    // cout << GWeight << ", " << VerWeight << endl;
    return GWeight * VerWeight * Factor;
  }

  // Sigma with LoopNum>=2
  G1.K = Var.LoopMom[1];
  G2.K = Var.LoopMom[2];
  G3.K = Var.LoopMom[1] + Var.LoopMom[2] - Var.LoopMom[0];
  G1.Evaluate();
  G2.Evaluate();
  G3.Evaluate();

  //   cout << "G: " << G1[0] << ", " << G2[0] << ", " << G3[0] << endl;

  double Weight = 0.0;

  for (auto &b : Bubble) {
    b.LVer.Evaluate(Var.LoopMom[0], G1.K, G3.K, G2.K, false);
    b.RVer.Evaluate(G2.K, G3.K, G1.K, Var.LoopMom[0], false);

    // cout << "left: " << b.LVer.Weight[DIR] << ", " << b.LVer.Weight[EX] <<
    // endl; cout << "righ: " << b.RVer.Weight[DIR] << ", " <<
    // b.RVer.Weight[EX]
    // << endl;

    for (auto &map : b.Map) {
      auto &LvW = b.LVer.Weight[map[0]];
      auto &RvW = b.RVer.Weight[0];
      double w = (LvW[DIR] * RvW[DIR] + LvW[EX] * RvW[EX]) * SPIN;
      //   cout << w << endl;
      w += LvW[DIR] * RvW[EX] + LvW[EX] * RvW[DIR];
      w *= G1[map[1]] * G2[map[2]] * G3[map[3]] * 0.5;
      //   cout << "G" << G1[map[1]] * G2[map[2]] * G3[map[3]] << endl;
      Weight += w;
    }
  }

  //   Weight =
  //       1.0 / pow((G1.K.squaredNorm() + 1), 2) / pow((G2.K.squaredNorm() +
  //       1), 2);

  return Weight * Factor * Factor;
}

string sigma::ToString() {
  if (Order <= 1)
    return fmt::format("Order {0} is empty", Order);
  string indent = "";
  string Info =
      indent + fmt::format("├Sigma: LoopNum: {0}, G1:{1}, G2: {2}, G3: {3}\n",
                           LoopNum(), G1.Size(), G2.Size(), G3.Size());
  Info += indent + fmt::format("└─\n");
  for (int p = 0; p < Bubble.size(); p++) {
    Info += indent + ". │\n";
    verPair &pp = Bubble[p];
    Info +=
        indent + fmt::format(". ├PAIR - LVerLoopNum: {0}\n", pp.LVer.LoopNum());
    Info += indent + fmt::format(". ├─Map:");

    ASSERT_ALLWAYS(
        pp.Map.size() == pp.LVer.Tpair.size() * pp.RVer.Tpair.size(),
        fmt::format("Map size {0} != LVer.Tpair size {1} *RVer.Tpair size {2}!",
                    pp.Map.size(), pp.LVer.Tpair.size(), pp.RVer.Tpair.size()));

    for (auto &m : pp.Map)
      Info += fmt::format("({0},{1},{2},{3}), ", m[0], m[1], m[2], m[3]);

    Info += "\n";
    // Info += "\n" + indent + ". │\n";
    Info += pp.LVer.ToString(indent + ". ");
    // Info +=
    //     indent +
    //     ".....................................................\n";

    // Info += indent + ". │\n";
    Info += pp.RVer.ToString(indent + ". ");
  }

  // Info += "=======================================================\n";
  // Info += indent +
  // "----------------------------------------------------\n";
  return Info;
}

bool sigma::Test() {
  //   return false;
  // Two-loop sigma
  if (Order != 2)
    return false;

  double Factor = 1.0 / pow(2.0 * PI, D);
  int ExtTauIdx = MaxTauNum - 1;
  momentum &ExtK = Var.LoopMom[0];
  momentum &K1 = Var.LoopMom[1];
  momentum &K2 = Var.LoopMom[2];
  momentum K3 = K1 + K2 - ExtK;
  double Tau = Var.Tau[ExtTauIdx] - Var.Tau[0];
  double G1 = Prop.Green(Tau, K1, UP, 0);
  double G2 = Prop.Green(Tau, K2, UP, 0);
  double G3 = Prop.Green(-Tau, K3, UP, 0);
  double VerDir = Prop.Interaction(K1 - ExtK, 0);
  double VerEx = -Prop.Interaction(K2 - ExtK, 0);
  double VerW = (VerDir * VerDir + VerEx * VerEx) * SPIN + 2.0 * VerDir * VerEx;
  double Weight = VerW * G1 * G2 * G3 * Factor * Factor * 0.5;

  //   cout << "Test G: " << G1 << ", " << G2 << ", " << G3 << endl;
  //   cout << "Test: " << VerDir << ", " << VerEx << endl;

  ASSERT_ALLWAYS(IsEqual(Weight, Evaluate()),
                 "Sigma weight error: " << Weight << " vs " << Evaluate());
  // cout << "G_test=" << G1 << ", " << G2 << ", " << G3 << endl;
  // cout << "Ver_test=" << VerWeightDir << ", " << VerWeightExLeft << endl;

  //   cout << "DIR=" << Weight1 + 2 * Weight2 << ", EX=" << -Weight3 << endl;
  //   return Weight1 + 2 * Weight2;
  // cout << "Pass" << endl;
  return true;
}