#ifndef weight_H
#define weight_H

#include "dse.h"
#include "utility/rng.h"
#include "utility/utility.h"
#include "vertex.h"
#include <vector>
extern parameter Para;
extern RandomFactory Random;

namespace diag {
using namespace std;

#define MAXMOMNUM get_power<2, MaxOrder + 1>::value * 4

struct variable {
  int CurrOrder;
  int CurrChannel; // 0: I, 1: T, 2: U, 3: S
  long int CurrVersion;
  diagram CurrDiagram;

  int CurrExtMomBin; // current bin of the external momentum
  double CurrTau;    // current external tau
  ver::weightMatrix CurrWeight;
  double CurrAbsWeight;

  array<momentum, MaxMomNum> LoopMom; // all momentum loop variables
  array<double, MaxTauNum> Tau;       // all tau variables

  // variational approach, interactions with counterterms
};

class weight {
public:
  variable Var; // The variable of the integral
  array<ver::weightMatrix, 4> ChanWeight;

  ver::fermi Fermi;         // propagator factory
  ver::verQTheta VerQTheta; // vertex factory

  void Initialization();

  double Evaluate(int LoopNum, diagram Diagram);
  double EvaluateUST(int LoopNum);

  void MeasureUST();

  void EvaluateChanVer4(int LoopNum, array<ver::weightMatrix, 4> ChanWeight);

  // initialization, read diagrams, then initialize variables

private:
  dse::verDiag VerDiag; // diagram factory
  // diagram for different order and channel
  dse::ver4 Ver4Root[MaxOrder];

  void Vertex4(dse::ver4 &Ver4, bool IsFast);

  void Ver0(dse::ver4 &Ver4);

  void ChanI(dse::ver4 &Ver4, bool IsFast);
  void ChanUST(dse::ver4 &Ver4, bool IsFast);
};

}; // namespace diag

#endif
