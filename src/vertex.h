#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"
#include <array>

double sum2(const momentum &);
double norm2(const momentum &);

namespace diag {

const int MAXSIGMABIN = 100000;

class fermi {
public:
  fermi();
  double Green(double Tau, const momentum &Momentum, spin Spin, int GType,
               double Scale = 0);

private:
  // beyond which the expense sigma function will be called
  double UpperBound, LowerBound;
  double DeltaK;
  double UpperBound2, LowerBound2; // lower upbound for better sigma
  double DeltaK2;
  double PhyGreen(double Tau, const momentum &Mom, int GType, double Scale = 0);
  double FockSigma(const momentum &Mom);
  double BuildFockSigma();
  double Fock(double k);
  // warning: this function only works for T=0 and 3D!!!!
  double GetSigma(double k);
  double Sigma[MAXSIGMABIN];
  double Sigma2[MAXSIGMABIN];
};

class verQTheta {
public:
  verQTheta();
  void Interaction(const array<momentum *, 4> &LegK, double Tau, int VerType,
                   double &WeightDir, double &WeightEx);

  void Measure(const momentum &InL, const momentum &InR, const int QIndex,
               int Order, double Tau, int Channel, double WeightFactor);
  void Update(double Ratio, int Order);
  void Save(bool Simple = false);
  void ClearStatis();
  void LoadWeight();
  void ResetIRScale(int IRScaleBin);
  double *ChanT;
  double *dChanT;
  double *intChanT;

  double *ChanU;
  double *dChanU;

  double *ChanS;
  double *dChanS;

  double *ChanI;
  double *dChanI;

  double &EffInterT(int Angle, int ExtQ);
  double &DiffInterT(int Order, int Angle, int ExtQ);

  double &EffInterU(int Angle, int ExtQ);
  double &DiffInterU(int Order, int Angle, int ExtQ);

  double &EffInterS(int Angle, int ExtQ);
  double &DiffInterS(int Order, int Angle, int ExtQ);

  double &EffInterI(int Angle, int ExtQ);
  double &DiffInterI(int Order, int Angle, int ExtQ);
  // double EffInteraction[ScaleBinSize + 1][AngBinSize][ExtMomBinSize];
  // double DiffInteraction[MaxOrder][ScaleBinSize +
  // 1][AngBinSize][ExtMomBinSize]; double IntInteraction[MaxOrder][ScaleBinSize
  // + 1][AngBinSize][ExtMomBinSize];

  // double TauBasis[TauBinSize][TauBasisNum];

  double Normalization;
  double PhyWeightT;
  double PhyWeightI;

  int QIndex;
  int AngleIndex;
  int OrderIndex;

  int QIndexI;
  int AngleIndexI;
  int OrderIndexI;
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

}; // namespace diag
#endif
