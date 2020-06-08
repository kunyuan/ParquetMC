#ifndef propagator_H
#define propagator_H
#include "global.h"
#include <Eigen/Dense>

using namespace std;

namespace diag {

using namespace Eigen;
typedef VectorXd weight1D;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> weight2D;

class propagator {
public:
  propagator();
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
  void TestF();
  void LoadGreen();
  double ExtrapF(double Tau,double K);
  


private:
  double _BareGreen(double Tau, const momentum &K, spin Spin, int GType = 0);

  double _Interp1D(double K, const weight1D &data);
  double _Interp2D(double K, double tau, const weight2D &data);
  vector<double> _f;
  vector<double> _taulist;
  vector<double> _extMom;
  weight1D _StaticSigma;
  weight2D _DeltaG;
};

double Angle3D(const momentum &K1, const momentum &K2);
void _TestAngle2D();

} // namespace diag

#endif
