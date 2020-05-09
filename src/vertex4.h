#ifndef vertex4_H
#define vertex4_H

#include "global.h"
#include "propagator.h"
#include "utility/utility.h"
#include <array>
#include <string>
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
  int AddTidxPair(const array<int, 2> &Tpair);           // Add T pair
  double operator[](int Tidx) { return _Weight[Tidx]; }; // get _Weight[Tidx]
  // evaluate all weight with different T pairs and a given K
  void Evaluate(bool IsAnomal = false);
  void Evaluate(const momentum &K, bool IsAnomal = false);
  int Size() { return _Tpair.size(); };

  vector<array<int, 2>> _Tpair;
  vector<double> _Weight;
};

struct bubble;

class vertex4 {
public:
  void Build(int Level, int Order, int LoopIdx, int InTIdx,
             const vector<channel> &Channel, int Side, bool InBox);

  vector<channel> Channel;      // list of channels except I
  vector<verWeight> ChanWeight; // the weight of each channel

  array<green, 6> G;           // G for 0, T, U, S, TC and UC
  vector<array<int, 4>> Tpair; // external T list
  vector<verWeight> Weight;    // size: equal to T.size()

  vector<bubble> _UST; // bubble diagrams

  int TauNum() { return Order + 1; }
  int LoopNum() { return Order; }
  bool InBox() { return _InBox; }
  void
  Evaluate(const momentum &KInL, const momentum &KOutL, const momentum &KInR,
           const momentum &KOutR,
           bool IsFast = false); // evaluate the weights with different Tidx

  string ToString(string indent = "");
  void Test() { _TestOneLoopGamma(); };

private:
  int Level;
  int Side; // right side vertex is always a full gamma4
  int Order;
  int Tidx;
  int LoopIdx;
  bool _InBox; // this Ver4 is a SUB-diagram within a box

  // vector<envelope> Envelope; // envelop diagrams

  void _BuildBare();
  void _BuildI(channel chan);
  // create UST bubble with a given channel and left vertex order
  bubble _BuildBubble(channel chan, int ol);

  // utility functions
  int _AddTidxPair(const array<int, 4> &T);

  void _EvalBare(const momentum &KInL, const momentum &KOutL,
                 const momentum &KInR, const momentum &KOutR);
  void _EvalUST(const momentum &KInL, const momentum &KOutL,
                const momentum &KInR, const momentum &KOutR,
                bool IsFast = false);
  void _EvalI(const momentum &KInL, const momentum &KOutL, const momentum &KInR,
              const momentum &KOutR, bool IsFast = false);

  void _TestOneLoopGamma();
  verWeight _GetWeight(int Order, vector<channel> Channel);
};

struct bubble {
  channel Channel;
  vertex4 LVer;
  vertex4 RVer;
  // map LVerIdx and RVerIdx to VerIdx
  // LVerIdx, RVerIdx, G0idx, Gidx, VerIdx
  vector<array<int, 5>> Map;
};

enum index { LVERT, RVERT, G0T, GXT, VERT };

} // namespace diag

#endif