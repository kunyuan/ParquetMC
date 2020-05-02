#include "vertex4.h"

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

int green::AddTidxPair(const Vector2i &T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < _Tidx.size(); ++i) {
    auto &t = _Tidx[i];
    if (t[OUT] == T[OUT] && t[IN] == T[IN])
      return i;
  }
  _Tidx.push_back(T);
  _Weight.push_back(0.0);
  return _Tidx.size() - 1;
}

void green::Evaluate(const momentum &K) {
  int Size = _Tidx.size();
  for (int i = 0; i < Size; ++i) {
    auto &T = _Tidx[i];
    _Weight[i] = Prop.Green(Var.Tau[T[OUT]] - Var.Tau[T[IN]], K, UP);
  }
}
