#ifndef diagram_H
#define diagram_H

#include "global.h"
#include "utility/utility.h"
#include "vertex.h"
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

extern parameter Para;

namespace diag {

enum caltype { BARE, RG, PARQUET, RENORMALIZED, VARIATIONAL };
enum channel { I = 0, T, U, S, SIGMA };

struct bubble;
struct envelope;

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

} // namespace diag

#endif