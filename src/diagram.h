#ifndef diagram_H
#define diagram_H

#include "global.h"
#include "utility/utility.h"
#include <array>
#include <string>
#include <vector>

extern parameter Para;

namespace diag {
using namespace std;
const size_t MaxBranchNum = 1 << (MaxOrder - 1); // 2**(MaxOrder-1)

// column-major two dimensional array
template <class T, size_t ROW, size_t COL>
using matrix = std::array<std::array<T, COL>, ROW>;

typedef std::array<double, MaxLoopNum>
    loop; // array to store the loop basis for a propagator or interaction line
typedef std::array<int, 2> tau; // array to store the tau basis (In and Out)
                                // for a propagator or interaction line

// A group can be diagrams with different orders,
// or diagrams with same order but have little sign cancellation
// store G and Ver indexes pointing to the corresponding pool
struct group {
  std::string Name;
  int ID;
  int HugenNum;           // Number of Hugenholz diagrams in each group
  int Order;              // diagram order of the group
  int Ver4Num;            // number of 4-vertex
  int GNum;               // number of G
  int LoopNum;            // dimension of loop basis
  int InternalLoopNum;    // dimension of internal loop basis
  int ExtLoopNum;         // dimension of external loop basis
  int ExtTransferLoopNum; // dimension of external loop basis
  int ExtLegLoopNum;      // dimension of external loop basis
  int TauNum;             // dimension of tau basis
  int ExtTauNum;          // dimension of external tau basis
  int InternalTauNum;     // dimension of internal tau basis
  double ReWeight;
  double Weight;
  double NewWeight;
  array<bool, MaxLoopNum> IsExtLoop;
  array<bool, MaxLoopNum> IsExtTransferLoop;
  array<bool, MaxLoopNum> IsExtLegLoop;
  array<bool, MaxLoopNum> IsLockedLoop;
  array<bool, MaxTauNum> IsExtTau;
  array<bool, MaxTauNum> IsLockedTau;
};

group ReadOneGroup(istream &);

void Test(group &);

}; // namespace diag

std::string ToString(const diag::group &);

#endif