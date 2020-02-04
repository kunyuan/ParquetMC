#ifndef markov_H
#define markov_H

#include "global.h"
#include "utility/rng.h"
#include "weight.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace mc {
using namespace std;
const int MCUpdates = 7;

typedef array<double, ExtMomBinSize> polar;

class markov {
public:
  markov();
  long long Counter;

  void PrintMCInfo();
  void PrintDeBugMCInfo();
  void AdjustGroupReWeight();

  // MC updates
  void ChangeTau();
  void ChangeMomentum();
  void ChangeOrder();
  void ChangeScale();
  void ChangeChannel();

  void Measure();
  void ClearStatis();
  void SaveToFile(bool Simple);
  void LoadFile();

  int DynamicTest();

  // MC variables
  diag::weight Weight;
  diag::variable &Var;

private:
  array<ver::weightMatrix, 2> NewWeight;
  double NewAbsWeight;

  int GetTauNum(int Order);
  int GetLoopNum(int Order);

  // MC updates

  double ShiftExtTransferK(const int &, int &);
  double ShiftExtLegK(const momentum &, momentum &);
  double ShiftK(const momentum &, momentum &);
  double ShiftTau(const double &, double &);

  double GetNewTau(double &);
  double GetNewK(momentum &);
  double RemoveOldTau(double &);
  double RemoveOldK(momentum &);

  // MC updates information
  std::string UpdatesName[MCUpdates];
  double Accepted[MCUpdates][MaxOrder + 1];
  double Proposed[MCUpdates][MaxOrder + 1];

  enum Updates {
    INCREASE_ORDER = 0,
    DECREASE_ORDER,
    CHANGE_GROUP,
    CHANGE_TAU,
    CHANGE_MOM,
    CHANGE_SCALE,
    CHANGE_CHANNEL,
    END
  };
  std::string _DetailBalanceStr(Updates op);
};
}; // namespace mc

#endif
