#ifndef weight_H
#define weight_H

#include "diagram.h"
#include "dse.h"
#include "observable.h"
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
  long int CurrVersion;

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
  weight() : SigData(Var.Tau) {}
  variable Var; // The variable of the integral
  array<ver::weightMatrix, 4> ChanWeight;

  ver::fermi Fermi;         // propagator factory
  ver::verQTheta VerQTheta; // vertex factory
  ver::sigData SigData;

  dse::ver4 Ver4Root[MaxOrder];
  dse::sigma Sigma[MaxOrder];
  dse::polar Polar[MaxOrder];

  void Initialization();

  double Evaluate(int LoopNum);

  void MeasureUST();
  void MeasureSigma();
  void MeasurePolar();

  // void EvaluateChanVer4(int LoopNum, array<ver::weightMatrix, 4> ChanWeight);

  // initialization, read diagrams, then initialize variables
  void Benchmark(int LoopNum, diagram Diagram, int Step);
  void Test(int LoopNum, diagram Diagram);

private:
  dse::verDiag VerDiag; // diagram factory
  // diagram for different order and channel

  double EvaluateGamma(int LoopNum);
  double EvaluateSigma(int LoopNum, bool IsFast = true);
  double EvaluatePolar(int LoopNum);

  void Vertex4(dse::ver4 &Ver4, bool IsFast);

  void Ver0(dse::ver4 &Ver4);

  void ChanI(dse::ver4 &Ver4, bool IsFast);
  void ChanUST(dse::ver4 &Ver4, bool IsFast);

  void EvaluateG(vector<dse::green> &G, const momentum &K);
  double EvaluateUST(int LoopNum);

  void _TestOneLoopGamma();
  void _TestTwoLoopGamma();
  ver::weightMatrix _GetWeight(int LoopNum, vector<dse::channel> Channel);
};

}; // namespace diag

#endif
