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

class green {
public:
  int AddTidxPair(const Vector2i &Tpair);                // Add T pair
  double operator[](int Tidx) { return _Weight[Tidx]; }; // get _Weight[Tidx]
  // evaluate all weight with different T pairs and a given K
  void Evaluate(const momentum &Mom);

private:
  vector<Vector2i> _Tidx;
  vector<double> _Weight;
};

} // namespace diag

#endif