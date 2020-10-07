#define FMT_HEADER_ONLY
#include "diagram.h"
#include "utility/fmt/format.h"
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

double legendre(double xi, int channel){
  if(channel==0) return 1;
  if(channel==1) return xi;
  if(channel==2) return 0.5*(3*xi*xi-1);
  else return 0;
}

double epsilon(double p){
  double E=abs(p*p-Para.Ef);
  // // if(E<=EPS){
  // //   return 4/Para.Beta;
  // // }
  // double result= E*(2+2*cosh(Para.Beta*E))/sinh(Para.Beta*E);
  // if(!std::isfinite(result)) return 4/Para.Beta;
  // return result;
  return 1.0;//4*E*Para.Beta+EPS;
}

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
  double result=0;
  double Factor = 1.0 / pow(2.0 * π, D);
  int main_channel=0;
  // normalization
  if (Order == 0)
    return 1.0;
  else if (Order == 1) {
    // bare interaction
    double Weight = Prop.Interaction(Var.LoopMom[0] - Var.LoopMom[1], -1);
    Weight *= Prop.F(1.0e-8, Var.LoopMom[1], UP, 0,main_channel);
    // cout << "1: " << Weight * Factor * 0.5 << endl;
    double xi=Var.LoopMom[0].dot(Var.LoopMom[1])/Var.LoopMom[0].norm()/Var.LoopMom[1].norm();
    result= Weight * Factor * legendre(xi,main_channel);
    // if(result<1e-300){
    //   string output=fmt::format("change_mom:{0:e}\t{1:e}",Prop.Interaction(Var.LoopMom[0] + Var.LoopMom[1], -1),Prop.F(1.0e-8, Var.LoopMom[1], UP, 0));
    //   cout<<output<<endl;
    //   throw std::invalid_argument("delta eval order 1");
    // }
    return result*epsilon(Var.LoopMom[0].norm())/epsilon(Var.LoopMom[1].norm());
  }

  // loop order >=2
  vertex4 &Ver4 = Vertex;
  F.Evaluate(Var.LoopMom[1], main_channel);
  // if (Var.CurrOrder == 2)
  //   cout << Var.LoopMom[1].norm() << endl;

  Vertex.Evaluate(Var.LoopMom[0], Var.LoopMom[1], -Var.LoopMom[0],
                  -Var.LoopMom[1], false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &fidx = Fidx[i];
    Weight += (Ver4.Weight[i][DIR] - Ver4.Weight[i][EX]) * F[fidx];
    //Weight += (Ver4.Weight[i][DIR]) * F[fidx];

  }
  // there is a symmetry factor -0.5
  // cout << "2: " << Weight * Factor * 0.5 << endl;
  double xi=Var.LoopMom[0].dot(Var.LoopMom[1])/Var.LoopMom[0].norm()/Var.LoopMom[1].norm();
  result= Weight * Factor * (0.5)* legendre(xi,main_channel);

  // if(result<1e-300){
  //     throw std::invalid_argument("delta Eval order 2");
  // }
  return result*epsilon(Var.LoopMom[0].norm())/epsilon(Var.LoopMom[1].norm());
}

double delta::Evaluate(int channel) {
  double Factor = 1.0 / pow(2.0 * PI, D);
  double result=0;
  // normalization
  if (Order == 0)
    return 1.0;
  else if (Order == 1) {
    // bare interaction
    double Weight = Prop.Interaction(Var.LoopMom[0] - Var.LoopMom[1], -1);
    Weight *= Prop.F(1.0e-8, Var.LoopMom[1], UP, 0,channel);
    // cout << "1: " << Weight * Factor * 0.5 << endl;
    double xi=Var.LoopMom[0].dot(Var.LoopMom[1])/Var.LoopMom[0].norm()/Var.LoopMom[1].norm();
    result= Weight * Factor * legendre(xi,channel);
    // if(!std::isfinite(result)){
    //   throw std::invalid_argument("delta_evaluate nan");
    // }
    return result*epsilon(Var.LoopMom[0].norm())/epsilon(Var.LoopMom[1].norm());
  }
  // loop order >=2
  vertex4 &Ver4 = Vertex;
  F.Evaluate(Var.LoopMom[1], channel);
  // if (Var.CurrOrder == 2)
  //   cout << Var.LoopMom[1].norm() << endl;

  Vertex.Evaluate(Var.LoopMom[0], Var.LoopMom[1], -Var.LoopMom[0],
                  -Var.LoopMom[1], false);

  int Size = Vertex.Tpair.size();
  double Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &fidx = Fidx[i];
    Weight += (Ver4.Weight[i][DIR] - Ver4.Weight[i][EX]) * F[fidx];
    //Weight += (Ver4.Weight[i][DIR]) * F[fidx];
  }
  // there is a symmetry factor -0.5
  // cout << "2: " << Weight * Factor * 0.5 << endl;
  double xi=Var.LoopMom[0].dot(Var.LoopMom[1])/Var.LoopMom[0].norm()/Var.LoopMom[1].norm();
  result= Weight * Factor * (0.5) * legendre(xi,channel);
  // if(!std::isfinite(result)){
  //   throw std::invalid_argument("delta_evaluate nan");
  // }
  return result*epsilon(Var.LoopMom[0].norm())/epsilon(Var.LoopMom[1].norm());
}

string delta::ToString() { return Vertex.ToString(""); }
