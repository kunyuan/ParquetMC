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

class weight {
public:
  weight() {}
  array<verWeight, 4> ChanWeight;

  ver::fermi Fermi;         // propagator factory
  ver::verQTheta VerQTheta; // vertex factory
  ver::sigData SigData;
  ver::polarData PolarData;
  ver::deltaData DeltaData;

  dse::ver4 Ver4Root[MaxOrder];
  dse::sigma Sigma[MaxOrder];
  dse::polar Polar[MaxOrder];
  dse::delta Delta[MaxOrder];

  void Initialization();

  double Evaluate(int LoopNum);

  void MeasureUST();
  void MeasureSigma();
  void MeasurePolar();
  void MeasureDelta();

  // void EvaluateChanVer4(int LoopNum, array<ver::weightMatrix, 4> ChanWeight);

  // initialization, read diagrams, then initialize variables
  void Benchmark(int LoopNum, int Step);
  void Test(int LoopNum);

private:
  dse::verDiag VerDiag; // diagram factory
  // diagram for different order and channel

  double EvaluateGamma(int LoopNum);
  double EvaluateSigma(int LoopNum, bool IsFast = true);
  double EvaluatePolar(int LoopNum);
  double EvaluateDelta(int LoopNum);

  void Vertex4(dse::ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
               const momentum &KInR, const momentum &KOutR, bool IsFast);

  void Ver0(dse::ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
            const momentum &KInR, const momentum &KOutR);

  void ChanI(dse::ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
             const momentum &KInR, const momentum &KOutR, bool IsFast);

  void ChanUST(dse::ver4 &Ver4, const momentum &KInL, const momentum &KOutL,
               const momentum &KInR, const momentum &KOutR, bool IsFast);

  void EvaluateG(vector<dse::green> &G, const momentum &K);
  double EvaluateUST(int LoopNum);

  void _TestOneLoopGamma();
  void _TestTwoLoopGamma();
  double _TestTwoLoopSigma();
  verWeight _GetWeight(int LoopNum, vector<dse::channel> Channel);
};

}; // namespace diag

#endif
