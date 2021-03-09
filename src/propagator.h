#ifndef propagator_H
#define propagator_H
#include "global.h"
#include <Eigen/Dense>
#include <map>

namespace diag {

using namespace Eigen;
typedef VectorXd weight1D;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> weight2D;

class propagator {
public:
  // propagator();
  void Initialize();
  double Green(double Tau, const momentum &K, spin Spin, int GType = 0);
  double F(double Tau, const momentum &K, spin Spin, int GType = 0);

  verWeight Interaction(const momentum &KInL, const momentum &KOutL,
                        const momentum &KInR, const momentum &KOutR,
                        double ExtQ = 0.0);

  verWeight InteractionTauBare(const momentum &KInL, const momentum &KOutL,
                               const momentum &KInR, const momentum &KOutR,
                               double inTL, double inTR, double ExtQ = 0.0);

  void InteractionTau(const momentum &KInL, const momentum &KOutL,
                      const momentum &KInR, const momentum &KOutR, double inTL,
                      double inTR, verWeight &w0, verWeight &w1, verWeight &w2,
                      double ExtQ = 0.0);

  // get the Direct part of the interaction
  double Interaction(const momentum &TranQ, int VerOrder = 0,
                     double ExtQ = 0.0);

  double CounterBubble(const momentum &K);

  void LoadF();
  void LoadGreen();
  void LoadGreenOrder();
  void SaveGreenOrder();

  void LoadInteraction();

private:
  std::map<int, weight1D> _StaticSigma;
  weight2D _DeltaG;
  weight2D _DeltaRs;
  weight2D _DeltaRa;
  std::map<int, weight2D> _deltaGOrder;
  grid::Tau _TauGridInterp;
  grid::FermiK _MomGridInterp;

  template <typename KGrid>
  double _Interp1D(const weight1D &data, const KGrid &kgrid, double K);
  template <typename KGrid>
  double _Interp2D(const weight2D &data, const KGrid &kgrid, double K,
                   double T);
  template <typename KGrid>
  double _Interp2D(const weight2D &data, const KGrid &kgrid,
                   const grid::Tau &tGrid, double K, double T);
};

double Angle3D(const momentum &K1, const momentum &K2);
void _TestAngle2D();

} // namespace diag

#endif