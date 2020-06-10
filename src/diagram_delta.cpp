#define FMT_HEADER_ONLY
#include "diagram.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void delta::Build(int order) {
  ASSERT_ALLWAYS(order >= 0, "Polar order must be larger than 0!");
  Order = order;

  // vertex is only needed for order>=2
  if (Order < 1)
    return;
  vector<channel> Chan = {I, T, U};
  // if the bare part of W is re-expaned, then TC and UC are also needed
  //   vector<channel> Chan = {I, T, U, TC, UC};
  Vertex.Build(0,         // level
               Order - 1, // loopNum
               2,         // loop index of the first internal K of the vertex
               0,         // tau index of the InTL leg
               Chan, RIGHT);
  for (auto &t : Vertex.Tpair) {
    int idx = F.AddTidxPair({t[OUTL], t[OUTR]});
    Fidx.push_back(idx);
  }
  // reset the last Tidx to MaxTanNum-1, which is discretized
  int ExtTidx = MaxTauNum - 1;
  ASSERT_ALLWAYS(ExtTidx > TauNum(), "MaxTauNum is too small!");

  for (auto &T : F._Tpair)
    for (auto &t : T)
      if (t == TauNum() - 1)
        t = ExtTidx;

  _ResetLastTidx(Vertex);
};

void delta::_ResetLastTidx(vertex4 &Vertex) {

  int LastTidx = TauNum() - 1;
  int ExtTidx = MaxTauNum - 1;

  for (auto &T : Vertex.Tpair)
    for (auto &t : T)
      if (t == LastTidx)
        t = ExtTidx;

  for (auto &g : Vertex.G)
    for (auto &T : g._Tpair)
      for (auto &t : T)
        if (t == LastTidx)
          t = ExtTidx;

  for (auto &bub : Vertex._UST) {
    _ResetLastTidx(bub.LVer);
    _ResetLastTidx(bub.RVer);
  }
  return;
}

double delta::Evaluate() {
  double Factor = 1.0 / pow(2.0 * π, D);
  // normalization
  if (Order == 0)
    return 1.0;
  else if (Order == 1) {
    // bare interaction
    double Weight = Prop.Interaction(Var.LoopMom[0] + Var.LoopMom[1], -1);
    Weight *= Prop.F(-1.0e-8, Var.LoopMom[1], UP, 0);
    // cout << "1: " << Weight * Factor * 0.5 << endl;
    return Weight * Factor;
  }

  // loop order >=2
  vertex4 &Ver4 = Vertex;
  F.Evaluate(Var.LoopMom[1], true);
  // if (Var.CurrOrder == 2)
  //   cout << Var.LoopMom[1].norm() << endl;

  Vertex.Evaluate(Var.LoopMom[0], Var.LoopMom[1], -Var.LoopMom[0],
                  -Var.LoopMom[1], false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &fidx = Fidx[i];
    Weight += (Ver4.Weight[i][DIR] - Ver4.Weight[i][EX]) * F[fidx];
  }
  // there is a symmetry factor -0.5
  // cout << "2: " << Weight * Factor * 0.5 << endl;
  return Weight * Factor * (-0.5);
}

string delta::ToString() { return Vertex.ToString(""); }