#ifndef dse_H
#define dse_H

#include "global.h"
#include "utility/utility.h"
#include "vertex.h"
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

extern parameter Para;

namespace dse {
using namespace std;

enum channel { I = 0, T, U, S, TC, UC };
const string ChanName[] = {"I", "T", "U", "S", "TC", "UC"};
const double SymFactor[6] = {1.0, -1.0, 1.0, -0.5, +1.0, -1.0};
const channel ChanMap[] = {I, T, U,
                           S, T, U}; // map channels to a more compact form

struct bubble;
struct envelope;

struct green {
  array<int, 2> T;
  double Weight;
};

struct ver4 {
  int Level;
  int Side; // right side vertex is always a full gamma4
  int TauNum;
  int InTL;
  int LoopNum;
  int Loopidx;
  bool InBox; // if this vertex has been the children of a boxed parent
              // vertex
  bool HasCT; // contain projected counter-term or not
  vector<channel> Channel; // list of channels except I
  // array<momentum, 4> ProjLegK; // projected external K derived from LegK

  /////////// bubble diagrams ////////////////////
  vector<bubble> Bubble; // different arrangement of LVer and RVer
  array<momentum, 4> K;  // momentum for internal K
  array<vector<green>, 6> G;
  // G lists for each channel, G0 is shared for all diagrams

  // vector<envelope> Envelope; // envelop diagrams and its counter diagram

  vector<array<int, 4>> T;  // external T list
  vector<verWeight> Weight; // size: equal to T.size()
};

struct indexMap {
  int LVerTidx; // LVer T index
  int RVerTidx; // RVer T index
  // map LVer T index and RVer T index to merged T index
  int Tidx;
  // LVer T and RVer T to Internal G0 and G
  int G0idx;
  int Gidx;
};

struct bubble {
  channel Channel;
  ver4 LVer;
  ver4 RVer;
  vector<indexMap> Map;
};

//////////////// Envelope diagrams /////////////////////////////

class g2Matrix {
public:
  g2Matrix() {
    InT = 0;
    OutT = 0;
  }
  g2Matrix(int _InT, int _OutT, momentum *k) {
    K = k;
    InT = _InT;
    OutT = _OutT;
  }

  int InT, OutT; // possible InT and OutT
  momentum *K;   // momentum on G
  double Weight;
};

struct mapT4 {
  int LDVerTidx;
  int RDVerTidx;
  int LUVerTidx;
  int RUVerTidx;
  // map LVer T index and RVer T index to merged T index
  array<int, 4> Tidx; // external T for four envelop diagrams
  // LVer T and RVer T to Internal T for G1 and G2
  array<array<int, 2>, 9> GT;
};

struct envelope {
  bool IsProjected;
  int InTL;
  array<momentum *, 4> LegK; // legK index
  array<ver4, 10> Ver;
  array<g2Matrix, 9> G;
  vector<mapT4> Map;
  array<double, 4> SymFactor;
};

////////////// Vertex Creation Class /////////////////////////////////
class verDiag {
public:
  ver4 Vertex(int Level, int LoopNum, int LoopIndex, int InTL,
              const vector<channel> &Channel, int Side, bool InBox);
  string ToString(const ver4 &Vertex, string indent = "");
  void ResetMomMap(ver4 &Ver4, const array<momentum *, 4> &LegK);
  // the LegK pointers of Ver4 must be reset after all memory allocation is
  // done.
  void Ver0(ver4 &Ver4);

private:
  void ChanI(ver4 &Ver4, const vector<channel> &Channel);
  void ChanUST(ver4 &Ver4, const vector<channel> &Channel);
};

bool verTest();

int AddToTList(vector<array<int, 4>> &TList, const array<int, 4> &T);

int AddToGList(vector<green> &GList, const array<int, 2> &T);

vector<indexMap> CreateIndexMap(ver4 &Ver4, const ver4 &LVer, const ver4 &RVer,
                                channel Chan);

} // namespace dse
#endif