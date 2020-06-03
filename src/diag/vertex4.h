#ifndef vertex4_H
#define vertex4_H

#include "global.h"
#include "propagator.h"
#include "utility/utility.h"
#include <array>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace tree {
using namespace std;

typedef array<int, 2> t2pair;
typedef array<int, 4> t4pair;
typedef Eigen::Array<int, MaxMomNum, 1> kvector;

enum channel { I = 0, T, U, S, TC, UC };
const string ChanName[] = {"I", "T", "U", "S", "TC", "UC", "v"};
const double SymFactor[7] = {1.0, -1.0, 1.0, -0.5, +1.0, -1.0, 1.0};
// map channels to a more compact form
const channel ChanMap[] = {I, T, U, S, T, U, I};

struct bubble;
struct envelope;

kvector KVector(int Loopidx);

class vertex4 {
public:
  void Build(int Level, const set<channel> &Channel,
             const array<kvector, 4> &LegK, int loopNum, int Loopidx,
             int InTIdx, int Side);

  set<channel> Channel;   // list of channels except I
  set<channel> ChannelCT; // list of counterterm channels except I

  array<kvector, 4> LegK;

  vector<bubble> UST; // bubble diagrams
  vector<t4pair> T4;  // external T list

  int TauNum() { return LoopNum + 1; }

  string ToString(string indent = "");
  void Test();

private:
  int Level;
  int Side; // right side vertex is always a full gamma4
  int LoopNum;
  int Tidx;
  int LoopIdx;

  // vector<envelope> Envelope; // envelop diagrams

  void _BuildBare(int Tidx);
  void _BuildI(channel chan);
  // create UST bubble with a given channel and left vertex order
  bubble _BuildBubble(channel chan, int ol);

  // utility functions
  int _AddTidxPair(const t4pair &T);
  string _ToString(const set<channel> &Chan);
};

struct bubble {
  channel Channel;
  vertex4 LVer;
  vertex4 RVer;
  kvector K0, Kx;
  vector<t2pair> G0, Gx; // t2pair of two G
  // map LVerIdx and RVerIdx to VerIdx
  // LVerIdx, RVerIdx, G0idx, Gidx, VerIdx
  vector<array<int, 5>> Map;
};

enum index { LVER, RVER, G0, GX, VER };

} // namespace tree

#endif