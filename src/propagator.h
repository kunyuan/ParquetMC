#ifndef propagator_H
#define propagator_H
#include "global.h"

namespace diag {
class propagator {
public:
  void Initialize();
  double Green(double Tau, const momentum &K, spin Spin, int GType = 0);

  verWeight Interaction(const momentum &KInL, const momentum &KOutL,
                        const momentum &KInR, const momentum &KOutR,
                        bool Boxed = false, double ExtQ = 0.0);

  // get the Direct part of the interaction
  double Interaction(const momentum &TranQ, int VerOrder = 0);

  double CounterBubble(const momentum &K);

private:
  double _BareGreen(double Tau, const momentum &K, spin Spin, int GType = 0);
};

} // namespace diag

#endif