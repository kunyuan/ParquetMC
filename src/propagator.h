#ifndef propagator_H
#define propagator_H
#include "global.h"
#include <Eigen/Dense>

namespace diag {

using namespace Eigen;
typedef VectorXd weight1D;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> weight2D;

class propagator {
public:
  void Initialize();
  double Green(double Tau, const momentum &K, spin Spin, int GType = 0);
  double F(double Tau, const momentum &K, spin Spin, int GType = 0);

  verWeight Interaction(const momentum &KInL, const momentum &KOutL,
                        const momentum &KInR, const momentum &KOutR,
                        double ExtQ = 0.0);

  // get the Direct part of the interaction
  double Interaction(const momentum &TranQ, int VerOrder = 0,
                     double ExtQ = 0.0);

  double CounterBubble(const momentum &K);

  void LoadF();
  void LoadGreen();

private:
  weight1D _StaticSigma;
  weight2D _DeltaG;
  template <typename KGrid>
  double _Interp1D(const weight1D &data, const KGrid &kgrid, double K);
  template <typename KGrid>
  double _Interp2D(const weight2D &data, const KGrid &kgrid, double K,
                   double T);
};

double Angle3D(const momentum &K1, const momentum &K2);
void _TestAngle2D();

} // namespace diag

#endif