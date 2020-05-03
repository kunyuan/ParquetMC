#ifndef vertex4_H
#define vertex4_H

#include "global.h"
#include "propagator.h"
#include "utility/utility.h"
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace diag {
using namespace std;

enum channel { I = 0, T, U, S, TC, UC };
const string ChanName[] = {"I", "T", "U", "S", "TC", "UC"};
const double SymFactor[6] = {1.0, -1.0, 1.0, -0.5, +1.0, -1.0};
// map channels to a more compact form
const channel ChanMap[] = {I, T, U, S, T, U};

struct bubble;
struct envelope;

// two-point vertex function
class green {
public:
  momentum K;
  int AddTidxPair(const Vector2i &Tpair);                // Add T pair
  double operator[](int Tidx) { return _Weight[Tidx]; }; // get _Weight[Tidx]
  // evaluate all weight with different T pairs and a given K
  void Evaluate();

private:
  vector<Vector2i> _Tidx;
  vector<double> _Weight;
};

struct bubble;

class vertex4 {
public:
  int Level;
  int Side; // right side vertex is always a full gamma4
  int Order;
  int InTidx;
  int Loopidx;
  bool InBox;              // boxed vertex or not
  vector<channel> Channel; // list of channels except I

  array<vector<green>, 6> G; // G for 0, T, U, S, TC and UC

  vector<array<int, 4>> T;  // external T list
  vector<verWeight> Weight; // size: equal to T.size()

  int TauNum() { return Order + 1; }
  int LoopNum() { return Order; }
  void Build(int Level, int LoopNum, int LoopIndex, int InTL,
             const vector<channel> &Channel, int Side, bool InBox);

private:
  vector<bubble> _UST; // bubble diagrams
  // vector<envelope> Envelope; // envelop diagrams
};

struct bubble {
  channel Channel;
  vertex4 LVer;
  vertex4 RVer;
  // map LVerIdx and RVerIdx to VerIdx
  // LVerIdx, RVerIdx, G0idx, Gidx, VerIdx
  vector<array<int, 5>> Map;
};

enum index { LVER, RVER, G0, GX, VER };

} // namespace diag

#endif