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

enum caltype { BARE, RG, PARQUET, RENORMALIZED, VARIATIONAL };
enum channel { I = 0, T, U, S };

struct bubble;
struct envelope;

struct ver4 {
  int ID;
  int LoopNum;
  int TauNum;
  int Side; // right side vertex is always a full gamma4

  bool IsFullVer4;
  bool RenormVer4;  // renormalize the current vertex
  bool RexpandBare; // reexpand the coupling in the left vertex

  vector<bubble> Bubble;     // bubble diagrams and its counter diagram
  vector<envelope> Envelope; // envelop diagrams and its counter diagram

  array<momentum *, 4> LegK;        // external legK index
  vector<array<int, 4>> T;          // external T list
  vector<ver::weightMatrix> Weight; // size: equal to T.size()
};

//////////////// Bubble diagrams /////////////////////////////

class gMatrix {
public:
  gMatrix() {
    _TauNum = 0;
    _InTL = 0;
  }
  gMatrix(int TauNum, int InTL, momentum *k) {
    _TauNum = TauNum;
    _InTL = InTL;
    K = k;
    _G.resize(TauNum * TauNum);
    for (auto &g : _G)
      g = 0.0;
  }
  double &operator()(int l, int r) {
    return _G[(l - _InTL) * _TauNum + r - _InTL];
  }

  double &operator()(const array<int, 2> &t) {
    return _G[(t[IN] - _InTL) * _TauNum + t[OUT] - _InTL];
  }

  momentum *K;

private:
  int _TauNum;
  int _InTL;
  vector<double> _G;
};

struct mapT2 {
  int LVerTidx; // LVer T index
  int RVerTidx; // RVer T index
  // map LVer T index and RVer T index to merged T index
  int Tidx; // three channels
  // LVer T and RVer T to Internal T for G1 and G2
  array<int, 2> G0T; // the shared G
  array<int, 2> GT;
};

struct pair {
  ver4 LVer;
  ver4 RVer;
  channel Channel;
  double SymFactor;
  vector<mapT2> Map;
};

struct bubble {
  int InTL;
  bool IsProjected;
  bool HasTU;
  bool HasS;
  vector<channel> Channel; // list of channels except I
  array<double, 4> ProjFactor;
  array<array<momentum *, 4>, 4> LegK; // legK index for different channel
  array<gMatrix, 4> G;
  vector<pair> Pair; // different Tau arrangement and channel
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
  ver4 Build(array<momentum, MaxMomNum> &loopmom, int LoopNum,
             vector<channel> Channel, caltype Type);
  string ToString(const ver4 &Vertex, string indent = "", int Level = 0);

private:
  int DiagNum = 0;
  int MomNum = MaxLoopNum;
  array<momentum, MaxMomNum> *LoopMom; // all momentum loop variables

  ver4 Vertex(array<momentum *, 4> LegK, int InTL, int LoopNum, int LoopIndex,
              vector<channel> Channel, int Side, bool RenormVer4,
              bool RexpandBare, bool IsFullVer4);

  ver4 Ver0(ver4 Ver4, int InTL);
  ver4 ChanI(ver4 Ver4, vector<channel> Channel, int InTL, int LoopNum,
             int LoopIndex, bool IsProjected = false);
  ver4 ChanUST(ver4 Ver4, vector<channel> Channel, int InTL, int LoopNum,
               int LoopIndex, bool IsProjected = false);
  momentum *NextMom();
};

bool verTest();

///////////////////   Self Energy  ///////////////////////////////////////

struct sigma {
  int ID;
  int LoopNum;
  int TauNum;
  bool RexpandBare; // reexpand the coupling in the left vertex
  int InTidx;       // external Tau index
  int OutTidx;
  momentum *LegK; // external legK index

  array<gMatrix, 3> G;
  ver4 RVer;

  double Weight; // size: equal to T.size()
};

////////////// Self Energy Creation Class ////////////////////////////////
class sigmaDiag {
public:
  sigma Build(array<momentum, MaxMomNum> &loopmom, int LoopNum, caltype Type);
  string ToString(string indent = "", int Level = 0);

private:
  int MomNum = MaxLoopNum;
  array<momentum, MaxMomNum> *LoopMom; // all momentum loop variables
  momentum *NextMom();
};

} // namespace dse
#endif