#ifndef propagator_H
#define propagator_H
#include "global.h"

namespace diag {
class propagator {
public:
  void Initialize();
  double Green(double Tau, const momentum &K, spin Spin, int GType = 0);
  double F(double TauIn, double TauOut, const momentum &K, spin Spin,
           int GType = 0);

  verWeight Interaction(const momentum &KInL, const momentum &KOutL,
                        const momentum &KInR, const momentum &KOutR,
                        bool Boxed = false, double ExtQ = 0.0);

  // get the Direct part of the interaction
  double Interaction(const momentum &TranQ, int VerOrder = 0);

  double CounterBubble(const momentum &K);

private:
  double _BareGreen(double Tau, const momentum &K, spin Spin, int GType = 0);
  void LoadF();
  void LoadGreen();
};

double Angle3D(const momentum &K1, const momentum &K2);
double Index2Angle(const int &Index, const int &AngleNum);
int Angle2Index(const double &Angle, const int &AngleNum);
void _TestAngleIndex();
void _TestAngle2D();

double Index2Mom(const int &Index);
int Mom2Index(const double &K);

double Index2Scale(const int &Index);
int Scale2Index(const double &Scale);

double Index2Tau(const int &Index);
int Tau2Index(const double &Tau);

} // namespace diag

#endif