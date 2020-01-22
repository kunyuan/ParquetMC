#ifndef weight_H
#define weight_H

#include "diagram.h"
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

  int CurrExtMomBin; // current bin of the external momentum
  double CurrTau;    // current external tau
  double CurrScale;  // Current (Reference) Scale: Index=1, ..., ScaleBinSize
  int CurrIRScaleBin;

  array<momentum, MaxMomNum> LoopMom; // all momentum loop variables
  array<double, MaxTauNum> Tau;       // all tau variables

  dse::weightMatrix CurrWeight;
  double CurrAbsWeight;
  // array<int, MaxMomNum> LoopSpin;     // all spin variables
};

class weight {
public:
  vector<group> Groups;
  variable Var; // The variable of the integral

  dse::weightMatrix Evaluate(int LoopNum, int Channel);

  // initialization, read diagrams, then initialize variables
  void ReadDiagrams();

  // initialization, read diagrams, then initialize variables
  // MC updates related operations
  // double ChangeTemperature(double NewBeta);
  void ChangeMom(group &, int Index);
  void ChangeTau(group &, int TauIndex);
  // two tau index on the two sides of interaction
  void ChangeGroup(group &, bool Forced = false);
  // recalculate the weights in one group
  double GetNewWeight(group &); // return the current weight
  void AcceptChange(group &);
  void RejectChange(group &);

  void Measure(double WeightFactor);
  void ClearStatis();
  void LoadWeight();
  void Save(bool Simple = false);

  // run test in MC updates
  int DynamicTest();

  // Test before MC
  int StaticTest();

  string DebugInfo(group &);

private:
  string _ErrMsg(string);

  void Initialization();

  template <typename... TS> string ERR(string format, TS... args);

  fermi Fermi;
  verQTheta VerQTheta;

  dse::verDiag VerDiag;
  // diagram for different order and channel
  dse::ver4 Ver4Root[MaxOrder][4];

  void Vertex4(dse::ver4 &Ver4);

  void Ver0(dse::ver4 &Ver4);

  void ChanI(dse::ver4 &Ver4);
  void ChanUST(dse::ver4 &Ver4);
};

}; // namespace diag

#endif
