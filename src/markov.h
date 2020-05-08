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
const int MCUpdates = 6;

class markov {
public:
  markov();

  diag::weight Weight; // Weight handler

  // information printer
  void PrintMCInfo();
  void PrintDeBugMCInfo();
  void AdjustGroupReWeight();

  // MC updates
  void ChangeTau();
  void ChangeMomentum();
  void ChangeExtMomentum();
  void ChangeExtTau();
  void ChangeOrder();
  void ChangeScale();

private:
  double NewAbsWeight;

  double ShiftExtTransferK(const int &, int &);
  double ShiftExtLegK(const int &, int &);
  double ShiftExtTau(const int &, int &);
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
    CHANGE_TAU,
    CHANGE_MOM,
    CHANGE_EXTMOM,
    CHANGE_EXTTAU,
    END
  };
  std::string _DetailBalanceStr(Updates op);
};
}; // namespace mc

#endif
