#include "vertex.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include "utility/utility.h"
#include <cmath>
#include <iostream>

using namespace ver;
using namespace std;

extern parameter Para;

verTensor::verTensor() {

  AngleIndex = ExtMomBinSize * 2;
  OrderIndex = ExtMomBinSize * AngBinSize * 2;

  _Interaction = new double[AngBinSize * ExtMomBinSize * 2];
  _Estimator = new double[MaxOrder * AngBinSize * ExtMomBinSize * 2];
}

verTensor::~verTensor() {
  delete[] _Interaction;
  delete[] _Estimator;
}

void verTensor::Initialize() {
  for (int inin = 0; inin < AngBinSize; ++inin)
    for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
      for (int tIndex = 0; tIndex < TauBinSize; ++tIndex) {
        Interaction(inin, qIndex, 0) = 0.0;
        Interaction(inin, qIndex, 1) = 0.0;
        for (int order = 0; order < MaxOrder; ++order) {
          Estimator(order, inin, qIndex, 0) = 0.0;
          Estimator(order, inin, qIndex, 1) = 0.0;
        }
      }
    }
}

double &verTensor::Interaction(int Angle, int ExtQ, int Dir) {
  return _Interaction[Angle * AngleIndex + ExtQ * 2 + Dir];
}

double &verTensor::Estimator(int Order, int Angle, int ExtQ, int Dir) {
  return _Estimator[Order * OrderIndex + Angle * AngleIndex + ExtQ * 2 + Dir];
}

// double norm2(const momentum &Mom) { return sqrt(sum2(Mom)); }
verQTheta::verQTheta() {

  //_TestAngle2D();
  //_TestAngleIndex();

  Normalization = 1.0e-10;
  // PhyWeight = Para.Kf * (1.0 - exp(-Para.MaxExtMom / Para.Kf)) * 4.0 * PI *
  // PI;
  // PhyWeight = (1.0 - exp(-Para.MaxExtMom / Para.Kf)) /
  //             (1.0 - exp(-Para.MaxExtMom / Para.Kf / ExtMomBinSize)) * 4.0 *
  //             PI * PI;
  // PhyWeight =
  //     1.0 / Para.Beta / Para.Beta * ExtMomBinSize * 2.0 * PI * Para.Kf * 4.0;

  PhyWeightT = AngBinSize;
  PhyWeightI = AngBinSize;
  // PhyWeight = 2.0 * PI / TauBinSize * 64;
  // PhyWeight = 1.0;

  for (auto &c : Chan)
    c.Initialize();
}

weightMatrix verQTheta::Interaction(const momentum &KInL, const momentum &KOutL,
                                    const momentum &KInR, const momentum &KOutR,
                                    double Tau, bool IsRenorm, bool Boxed,
                                    double ExtQ) {

  weightMatrix Weight;

  double kDiQ = (KInL - KOutL).norm();
  Weight[DIR] =
      -8.0 * PI * Para.Charge2 / (kDiQ * kDiQ + Para.Mass2 + Para.Lambda);

  if (DiagType == SIGMA && abs(kDiQ) < 1.0e-10)
    Weight[DIR] = 0.0;

  // check irreducibility
  if (DiagType == POLAR && abs(kDiQ - ExtQ) < 1.0e-10)
    Weight[DIR] = 0.0;

  if (!Boxed) {
    double kExQ = (KInL - KOutR).norm();
    Weight[EX] =
        8.0 * PI * Para.Charge2 / (kExQ * kExQ + Para.Mass2 + Para.Lambda);

    if (DiagType == SIGMA && abs(kExQ) < 1.0e-10)
      Weight[EX] = 0.0;

    // check irreducibility
    if (DiagType == POLAR && abs(kExQ - ExtQ) < 1.0e-10)
      Weight[EX] = 0.0;

  } else
    Weight[EX] = 0.0;

  // cout << "Ver0: " << Weight[DIR] << ", " << Weight[EX] << endl;
  // cout << "extnal: " << ExtQ << ", " << kDiQ << endl;
  return Weight;
}

double verQTheta::Interaction(const momentum &TranQ) {
  double kQ = TranQ.norm();
  return -8.0 * PI * Para.Charge2 / (kQ * kQ + Para.Mass2 + Para.Lambda);
  // return pow(-8.0 * PI * Para.Charge2 / (kQ * kQ + Para.Mass2 + Para.Lambda),
  //            2) *
  //        Para.Lambda / 8.0 / PI;
}

void verQTheta::Measure(const momentum &InL, const momentum &InR,
                        const int QIndex, int Order,
                        const array<ver::weightMatrix, 4> &Weight,
                        double Factor) {
  // cout << Order << ", " << DiagNum << endl;
  if (Order == 0) {
    Normalization += 1.0 * Factor;
    // Normalization += WeightFactor;
  } else {
    // double Factor = 1.0 / pow(2.0 * PI, 2 * Order);
    double CosAng = Angle3D(InL, InR);
    // cout << InL.norm() << ", " << InR.norm() << endl;
    // cout << "(" << InR[0] << ", " << InR[1] << ", " << InR[2] << ")" << endl;

    int AngleIndex = Angle2Index(CosAng, AngBinSize);
    for (int chan = 0; chan < 4; ++chan) {
      Chan[chan].Estimator(Order, AngleIndex, QIndex, DIR) +=
          Weight[chan][DIR] * Factor;
      Chan[chan].Estimator(0, AngleIndex, QIndex, DIR) +=
          Weight[chan][DIR] * Factor;

      Chan[chan].Estimator(Order, AngleIndex, QIndex, EX) +=
          Weight[chan][EX] * Factor;
      Chan[chan].Estimator(0, AngleIndex, QIndex, EX) +=
          Weight[chan][EX] * Factor;
    }
  }
  return;
}

void verQTheta::Save(bool Simple) {

  for (int chan = 0; chan < 4; chan++) {
    for (int order = 0; order <= Para.Order; order++) {
      if (Simple == true)
        if (order != 0)
          continue;
      string FileName =
          fmt::format("vertex{0}_{1}_pid{2}.dat", order, chan, Para.PID);
      ofstream VerFile;
      VerFile.open(FileName, ios::out | ios::trunc);

      if (VerFile.is_open()) {

        VerFile << fmt::sprintf("#PID:%d, rs:%.3f, Beta: %.3f, Step: %d\n",
                                Para.PID, Para.Rs, Para.Beta, Para.Counter);

        VerFile << "# Norm: " << Normalization << endl;

        VerFile << "# AngleTable: ";
        for (int angle = 0; angle < AngBinSize; ++angle)
          VerFile << Para.AngleTable[angle] << " ";

        VerFile << endl;
        VerFile << "# ExtMomBinTable: ";
        for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
          VerFile << Para.ExtMomTable[qindex][0] << " ";
        VerFile << endl;

        for (int angle = 0; angle < AngBinSize; ++angle)
          for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
            for (int dir = 0; dir < 2; ++dir)
              VerFile << Chan[chan].Estimator(order, angle, qindex, dir) *
                             PhyWeightT
                      << "  ";
        VerFile.close();
      } else {
        LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
      }
    }
  }
}

void verQTheta::ClearStatis() {
  Normalization = 1.0e-10;
  for (int inin = 0; inin < AngBinSize; ++inin)
    for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
      double k = Index2Mom(qIndex);
      for (int order = 0; order < MaxOrder; ++order)
        for (int dir = 0; dir < 2; ++dir)
          for (auto &c : Chan)
            c.Estimator(order, inin, qIndex, dir) = 0.0;
    }
}

void verQTheta::LoadWeight() {
  try {
    for (int chan = 0; chan < 4; chan++) {
      string FileName = fmt::format("../weight{0}.data", chan);
      ifstream VerFile;
      VerFile.open(FileName, ios::in);
      if (VerFile.is_open()) {
        for (int angle = 0; angle < AngBinSize; ++angle)
          for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
            for (int dir = 0; dir < 2; ++dir)
              VerFile >> Chan[chan].Interaction(angle, qindex, dir);
        VerFile.close();
      }
    }
  } catch (int e) {
    LOG_INFO("Can not load weight file!");
  }
}

fermi::fermi() {
  UpperBound = 5.0 * Para.Ef;
  LowerBound = 0.0;
  DeltaK = UpperBound / MAXSIGMABIN;
  UpperBound2 = 1.2 * Para.Ef;
  LowerBound2 = 0.8 * Para.Ef;
  DeltaK2 = UpperBound2 / MAXSIGMABIN;

  // Sigma = MatrixXd::Zero(SigmaMomBinSize, SigmaTauBinSize);
}
double fermi::Fock(double k) {
  // warning: this function only works for T=0!!!!
  double l = sqrt(Para.Mass2 + Para.Lambda);
  double kF = Para.Kf;
  double fock = 1.0 + l / kF * atan((k - kF) / l);
  fock -= l / kF * atan((k + kF) / l);
  fock -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
          log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
  fock *= (-2.0 * kF) / PI;

  double shift = 1.0 - l / kF * atan(2.0 * kF / l);
  shift -= l * l / 4.0 / kF / kF * log(l * l / (l * l + 4.0 * kF * kF));
  shift *= (-2.0 * kF) / PI;

  return fock - shift;
}

double fermi::BuildFockSigma() {
  ASSERT_ALLWAYS(D == 3, "The Fock self energy is for 3D!");
  double fock, k;
  for (int i = 0; i < MAXSIGMABIN; ++i) {
    // k: (0^+, UpBound^-)
    // i=0 ==> k==0.5*DeltaK
    // i=MAXSIGMABIN-1 ==> k==(MAXSIGMABIN-0.5)*DeltaK
    k = (i + 0.5) * DeltaK + LowerBound;
    Sigma[i] = Fock(k);
    if (i > 0 && k <= LowerBound2 && k >= UpperBound2) {
      ASSERT_ALLWAYS(
          Equal(Sigma[i - 1], Sigma[i], 5.0e-5),
          fmt::format("Fock are not accurate enough! At k={0}: {1} vs {2}\n", k,
                      Sigma[i - 1], Sigma[i]));
    }
    // cout << k << " : " << Sigma[i] << " vs " << Fock(k) << endl;
  }

  for (int i = 0; i < MAXSIGMABIN; ++i) {
    // k: (0^+, UpBound^-)
    // i=0 ==> k==0.5*DeltaK
    // i=MAXSIGMABIN-1 ==> k==(MAXSIGMABIN-0.5)*DeltaK
    k = (i + 0.5) * DeltaK2 + LowerBound2;
    Sigma2[i] = Fock(k);
    if (i > 0) {
      ASSERT_ALLWAYS(Equal(Sigma2[i - 1], Sigma2[i], 5.0e-5),
                     fmt::format("The 2rd level Fock are not accurate enough!"
                                 "level! At k={0}: {1} vs {2}\n",
                                 k, Sigma2[i - 1], Sigma2[i]));
    }
    // cout << k << " : " << Sigma[i] << " vs " << Fock(k) << endl;
  }
};

double fermi::FockSigma(const momentum &Mom) {
  double k = Mom.norm(); // bare propagator
  double fock;
  if (k >= LowerBound2 && k < UpperBound2) {
    int i = (k - LowerBound2) / DeltaK2;
    fock = Sigma2[i];
  } else if ((k >= LowerBound && k < LowerBound2) ||
             (k >= UpperBound2 && k < UpperBound)) {
    int i = (k - LowerBound) / DeltaK;
    fock = Sigma[i];
  } else {
    fock = Fock(k);
  }
  // ASSERT_ALLWAYS(
  //     Equal(fock, Fock(k), 5.0e-5),
  //     fmt::format("Fock are not accurate enough! At k={0}: {1} vs {2}\n",
  //     k,
  //                 fock, Fock(k)));
  return fock + k * k;
}

double fermi::PhyGreen(double Tau, const momentum &Mom, int GType,
                       double Scale) {
  // return 1.0;
  // if tau is exactly zero, set tau=0^-
  double green, Ek, kk, k;
  if (Tau == 0.0) {
    Tau = -1.0e-10;
    // return EPS;
  }
  // equal time green's function
  if (GType == 1)
    Tau = -1.0e-10;

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }

  k = Mom.norm();
  Ek = k * k; // bare propagator
  // Ek += Fock(k);
  // Ek = FockSigma(Mom);

  //// enforce an UV cutoff for the Green's function ////////
  // if(Ek>8.0*EF) then
  //   PhyGreen=0.0
  //   return
  // endif

  double x = Para.Beta * (Ek - Para.Mu) / 2.0;
  double y = 2.0 * Tau / Para.Beta - 1.0;
  if (x > 100.0)
    green = exp(-x * (y + 1.0));
  else if (x < -100.0)
    green = exp(x * (1.0 - y));
  else
    green = exp(-x * y) / (2.0 * cosh(x));

  green *= s;

  // if (std::isnan(green))
  //   ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
  //                 << ", Ek=" << Ek << ", Green=" << green << ", Mom"
  //                 << ToString(Mom));
  return green;
}

double fermi::Green(double Tau, const momentum &Mom, spin Spin, int GType,
                    double Scale) {
  double green;
  if (GType >= 0) {
    green = PhyGreen(Tau, Mom, GType, Scale);
  } else if (GType == -1) {
    // green = PhyGreen(Tau, Mom, Scale);
    green = 1.0;

  } else if (GType == -2) {
    // Lower Scale Green's function
    Scale -= 1;
    green = PhyGreen(Tau, Mom, GType, Scale);
  } else {
    ABORT("GType " << GType << " has not yet been implemented!");
    // return FakeGreen(Tau, Mom);
  }
  return green;
}

double ver::Index2Mom(const int &Index) {
  return (Index + 0.5) / ExtMomBinSize * Para.MaxExtMom;
};

int ver::Mom2Index(const double &K) {
  return int(K / Para.MaxExtMom * ExtMomBinSize);
};

double ver::Angle3D(const momentum &K1, const momentum &K2) {
  // Returns the angle in radians between vectors 'K1' and 'K2'
  double dotp = K1.dot(K2);
  double Angle2D = dotp / K1.norm() / K2.norm();
  return Angle2D;
}

double ver::Index2Angle(const int &Index, const int &AngleNum) {
  // Map index [0...AngleNum-1] to the theta range [0.0, 2*pi)
  return (Index + 0.5) * 2.0 / AngleNum - 1.0;
}

int ver::Angle2Index(const double &Angle, const int &AngleNum) {
  // Map theta range  [0.0, 2*pi) to index [0...AngleNum-1]
  // double dAngle = 2.0 * PI / AngleNum;
  // if (Angle >= 2.0 * PI - dAngle / 2.0 || Angle < dAngle / 2.0)
  //   return 0;
  // else
  //   return int(Angle / dAngle + 0.5);
  // if (Angle > 1.0 - EPS)
  //   return AngleNum - 1;
  // else {
  double dAngle = 2.0 / AngleNum;
  return int((Angle + 1.0) / dAngle);
  // }
}

int ver::Tau2Index(const double &Tau) {
  return int((Tau / Para.Beta) * TauBinSize);
}

double ver::Index2Tau(const int &Index) {
  return (Index + 0.5) * Para.Beta / TauBinSize;
}

void ver::_TestAngle2D() {
  // Test Angle functions
  momentum K1 = {1.0, 0.0};
  momentum K2 = {1.0, 0.0};

  ASSERT_ALLWAYS(
      abs(Angle3D(K1, K2) - 1.0) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not zero! It is {:.13f}",
                  Angle3D(K1, K2)));

  K1 = {1.0, 0.0};
  K2 = {-1.0, 0.0};
  ASSERT_ALLWAYS(
      abs(Angle3D(K1, K2) - (-1.0)) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not Pi! Instead, it is {:.13f}",
                  Angle3D(K1, K2)));

  K1 = {1.0, 0.0};
  K2 = {1.0, -EPS};
  ASSERT_ALLWAYS(
      abs(Angle3D(K1, K2) - 1.0) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not 2.0*Pi! It is {:.13f}",
                  Angle3D(K1, K2)));
}

void ver::_TestAngleIndex() {
  // Test Angle functions
  int AngleNum = 8;
  cout << Index2Angle(0, AngleNum) << endl;
  ASSERT_ALLWAYS(abs(Index2Angle(0, AngleNum) - (-1.0 + 1.0 / AngleNum)) <
                     1.0e-10,
                 "Angle for index 0 should be -1^+!");

  ASSERT_ALLWAYS(abs(Index2Angle(AngleNum - 1, AngleNum) -
                     (1.0 - 1.0 / AngleNum)) < 1.0e-10,
                 "Angle for index AngleNum should be 1.0^-!");

  ASSERT_ALLWAYS(Angle2Index(1.0 - EPS, AngleNum) == AngleNum - 1,
                 "cos(angle)=-1 should should be last element!");
  // ASSERT_ALLWAYS(
  //     Angle2Index(2.0 * PI * (1.0 - 0.5 / AngleNum) + EPS, AngleNum) == 0,
  //     "Angle 2*pi-pi/AngleNum should have index 1!");

  // ASSERT_ALLWAYS(Angle2Index(2.0 * PI * (1.0 - 0.5 / AngleNum) - EPS,
  //                            AngleNum) == AngleNum - 1,
  //                "Angle 2*pi-pi/AngleNum-0^+ should have index AngleNum!");
}
