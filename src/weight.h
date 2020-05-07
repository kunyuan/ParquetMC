#ifndef weight_H
#define weight_H

#include "diagram.h"
#include "observable.h"
#include "utility/utility.h"
#include "vertex4.h"
#include <vector>
// extern parameter Para;
// extern RandomFactory Random;

namespace diag {
using namespace std;

class weight {
public:
  weight() {}

  void Initialization();
  double Evaluate(int LoopNum);
  void Measure();
  void SaveToFile();
  void LoadFile(){};

  void Benchmark(int LoopNum, int Step);
  void Test();

  diag::vertex4 Gamma[MaxOrder];
  obs::ver4Obs GammaObs;

  diag::polar Polar[MaxOrder];
  diag::sigma Sigma[MaxOrder];
  obs::oneBodyObs OneBodyObs;
  // sigma, polar and delta can use the same observable

private:
};

}; // namespace diag

#endif
