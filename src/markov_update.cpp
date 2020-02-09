//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#include "markov.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

extern parameter Para;
extern RandomFactory Random;

using namespace mc;
using namespace diag;
using namespace std;

int markov::GetTauNum(int Order, int Channel) {
  if (Channel == dse::SIGMA)
    return Order + 2;
  else
    return Order + 1;
}

int markov::GetLoopNum(int Order) { return Order + 3; }

void markov::ChangeOrder() {
  Updates Name;
  double Prop = 1.0;
  int NewOrder;

  if (Random.urn() < 0.5) {
    // increase order
    Name = INCREASE_ORDER;
    if (Var.CurrOrder == Para.Order)
      return;
    NewOrder = Var.CurrOrder + 1;
    static momentum NewMom;
    double NewTau;
    // Generate New Tau
    Prop = GetNewTau(NewTau);
    int NewTauIndex = GetTauNum(Var.CurrOrder, Var.CurrChannel);
    Var.Tau[NewTauIndex] = NewTau;
    // Generate New Mom
    Prop *= GetNewK(NewMom);
    Var.LoopMom[GetLoopNum(Var.CurrOrder)] = NewMom;

    // if the current order is zero, then set channel of order 1 at T
    if (Var.CurrOrder == 0)
      Var.CurrChannel = dse::T;

  } else {
    // decrese order
    if (Var.CurrOrder == 0)
      return;

    // if the current order is one, then decrease order is possible only for T
    if (Var.CurrOrder == 1 && Var.CurrChannel != dse::T)
      return;

    Name = DECREASE_ORDER;
    NewOrder = Var.CurrOrder - 1;

    // Remove OldTau
    int TauToRemove = GetTauNum(Var.CurrOrder, Var.CurrChannel) - 1;
    Prop = RemoveOldTau(Var.Tau[TauToRemove]);
    // Remove OldMom
    int LoopToRemove = GetLoopNum(Var.CurrOrder) - 1;
    Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);
  }
  Proposed[Name][Var.CurrOrder] += 1;

  // Weight.ChangeGroup(NewGroup);
  NewWeight = Weight.Evaluate(NewOrder, Var.CurrChannel);
  NewAbsWeight = NewWeight.Abs();
  double R = Prop * NewAbsWeight * Para.ReWeight[NewOrder] / Var.CurrAbsWeight /
             Para.ReWeight[Var.CurrOrder];

  // if (Name == INCREASE_ORDER) {
  //   cout << "time:" << Var.Tau[2] << ", " << Var.Tau[3] << endl;
  // cout << NewWeight << endl;
  // }
  // if (NewGroup.Order == 2)
  //   cout << NewWeight << endl;

  if (Random.urn() < R) {
    Accepted[Name][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrOrder = NewOrder;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
  }
  return;
};

void markov::ChangeTau() {
  int TauIndex = Random.irn(0, GetTauNum(Var.CurrOrder, Var.CurrChannel) - 1);

  Proposed[CHANGE_TAU][Var.CurrOrder]++;

  double CurrTau = Var.Tau[TauIndex];
  double NewTau;
  double Prop = ShiftTau(CurrTau, NewTau);

  Var.Tau[TauIndex] = NewTau;

  NewWeight = Weight.Evaluate(Var.CurrOrder, Var.CurrChannel);
  NewAbsWeight = NewWeight.Abs();
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_TAU][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
  } else {
    // retore the old Tau if the update is rejected
    // if TauIndex is external, then its partner can be different
    Var.Tau[TauIndex] = CurrTau;
  }
};

void markov::ChangeMomentum() {
  int LoopIndex = Random.irn(0, GetLoopNum(Var.CurrOrder) - 1);
  // int LoopIndex = int(Random.urn() * (Var.CurrGroup->Order + 3));

  // InL momentum is locked
  if (LoopIndex == 1)
    return;

  Proposed[CHANGE_MOM][Var.CurrOrder]++;

  double Prop;
  int NewExtMomBin;
  static momentum CurrMom;

  CurrMom = Var.LoopMom[LoopIndex];

  if (LoopIndex == 0) {
    // transfer momentum
    Prop = ShiftExtTransferK(Var.CurrExtMomBin, NewExtMomBin);
    Var.LoopMom[LoopIndex] = Para.ExtMomTable[NewExtMomBin];
    if (Var.LoopMom[LoopIndex].norm() > Para.MaxExtMom) {
      Var.LoopMom[LoopIndex] = CurrMom;
      return;
    }
  } else if (LoopIndex == 2) {
    // InR momentum
    Prop = ShiftExtLegK(CurrMom, Var.LoopMom[LoopIndex]);
  } else {
    Prop = ShiftK(CurrMom, Var.LoopMom[LoopIndex]);
  }

  NewWeight = Weight.Evaluate(Var.CurrOrder, Var.CurrChannel);
  NewAbsWeight = NewWeight.Abs();
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_MOM][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
    if (LoopIndex == 0)
      Var.CurrExtMomBin = NewExtMomBin;
  } else {
    Var.LoopMom[LoopIndex] = CurrMom;
  }
};

void markov::ChangeChannel() {
  if (Var.CurrOrder == 0)
    return;
  double Prop = 1.0;
  int NewChannel = int(Random.urn() * 4);
  // Var.CurrChannel = dse::T;
  // if (Var.CurrChannel == dse::U) {
  //   Var.CurrChannel = OldChannel;
  //   return;
  // }
  Updates Name;

  if (Var.CurrChannel != dse::SIGMA && NewChannel == dse::SIGMA)
    Name = VER2SIGMA;
  else if (Var.CurrChannel == dse::SIGMA && NewChannel != dse::SIGMA)
    Name = SIGMA2VER;
  else
    Name = VER2VER;

  Proposed[Name][Var.CurrOrder] += 1;

  NewWeight = Weight.Evaluate(Var.CurrOrder, NewChannel);
  NewAbsWeight = NewWeight.Abs();
  double R = Prop * NewAbsWeight * Para.ReWeightChan[NewChannel] /
             Var.CurrAbsWeight / Para.ReWeightChan[Var.CurrChannel];
  if (Random.urn() < R) {
    Accepted[Name][Var.CurrOrder]++;
    Var.CurrVersion++;
    Var.CurrWeight = NewWeight;
    Var.CurrAbsWeight = NewAbsWeight;
    Var.CurrChannel = NewChannel;
  }
  return;
}

double markov::GetNewTau(double &NewTau) {
  double Step = 1.0;
  NewTau = Random.urn() * Para.Beta;
  // NewTau2 = (Random.urn() - 0.5) * Step + NewTau1;
  // if (NewTau2 < 0.0)
  //   NewTau2 += Para.Beta;
  // if (NewTau2 > Para.Beta)
  //   NewTau2 -= Para.Beta;
  // return Para.Beta * Step;
  return Para.Beta;
}

double markov::RemoveOldTau(double &OldTau) {
  // double Step = 1.0;
  // if (abs(OldTau2 - OldTau1) > Step / 2.0)
  //   return 0.0;
  // else
  //   return 1.0 / Para.Beta / Step;
  return 1.0 / Para.Beta;
}

double markov::GetNewK(momentum &NewMom) {
  //====== The hard Way ======================//
  // double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  // if (dK > Para.Kf / 2)
  // dK = Para.Kf / 2; // to avoid dK>Kf situation
  // double dK = 1.0 * Var.CurrScale;
  // double x = Random.urn();
  // double KAmp = dK / 2.0 + dK * x;
  // if (x < 0)
  //   KAmp = Para.Kf + KAmp;
  // else
  //   KAmp = Para.Kf - KAmp;
  double dK = Para.Kf / 2.0;
  double KAmp = Para.Kf + (Random.urn() - 0.5) * 2.0 * dK;
  if (KAmp <= 0.0) {
    return 0.0;
  }
  // Kf-dK<KAmp<Kf+dK
  double Phi = 2.0 * PI * Random.urn();
  if (D == 3) {
    double Theta = PI * Random.urn();
    if (Theta == 0.0)
      return 0.0;
    double K_XY = KAmp * sin(Theta);
    NewMom[0] = K_XY * cos(Phi);
    NewMom[1] = K_XY * sin(Phi);
    NewMom[D - 1] = KAmp * cos(Theta);
    return 2.0 * dK                    // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI                  // prop density of Phi
           * PI                        // prop density of Theta
           * sin(Theta) * KAmp * KAmp; // Jacobian
  } else if (D == 2) {
    NewMom[0] = KAmp * cos(Phi);
    NewMom[1] = KAmp * sin(Phi);
    return 2.0 * dK   // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI // prop density of Phi
           * KAmp;    // Jacobian
  }

  //===== The simple way  =======================//
  // for (int i = 0; i < D; i++)
  //   NewMom[i] = Para.Kf * (Random.urn() - 0.5) * 2.0;
  // return pow(2.0 * Para.Kf, D);

  //============================================//
};

double markov::RemoveOldK(momentum &OldMom) {
  //====== The hard Way ======================//
  // double dK = 1.0 * Var.CurrScale;
  // double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  // if (dK > Para.Kf / 2)
  //   dK = Para.Kf / 2; // to avoid dK>Kf situation
  double dK = Para.Kf / 2.0;
  double KAmp = OldMom.norm();
  if (KAmp < Para.Kf - dK || KAmp > Para.Kf + dK)
    // if (KAmp < Para.Kf - 1.5 * dK ||
    //     (KAmp > Para.Kf - 0.5 * dK && KAmp < Para.Kf + 0.5 * dK) ||
    //     KAmp > Para.Kf + 1.5 * dK)
    // Kf-dK<KAmp<Kf+dK
    return 0.0;
  if (D == 3) {
    auto SinTheta = sqrt(OldMom[0] * OldMom[0] + OldMom[1] * OldMom[1]) / KAmp;
    if (SinTheta < EPS)
      return 0.0;
    return 1.0 / (2.0 * dK * 2.0 * PI * PI * SinTheta * KAmp * KAmp);
  } else if (D == 2) {
    return 1.0 / (2.0 * dK * 2.0 * PI * KAmp);
  }

  //===== The simple way  =======================//
  // for (int i = 0; i < D; i++)
  //   if (fabs(OldMom[i] > Para.Kf))
  //     return 0.0;
  // return 1.0 / pow(2.0 * Para.Kf, D);
  //============================================//
}

double markov::ShiftK(const momentum &OldMom, momentum &NewMom) {
  double x = Random.urn();
  double Prop;
  if (x < 1.0 / 3) {
    // COPYFROMTO(OldMom, NewMom);
    NewMom = OldMom;
    int dir = Random.irn(0, D - 1);
    double STEP = Para.Beta > 1.0 ? Para.Kf / Para.Beta * 3.0 : Para.Kf;
    NewMom[dir] += STEP * (Random.urn() - 0.5);
    Prop = 1.0;
  } else if (x < 2.0 / 3) {
    double k = OldMom.norm();
    if (k < 1.0e-9) {
      Prop = 0.0;
      NewMom = OldMom;
    } else {
      const double Lambda = 1.5;
      double knew = k / Lambda + Random.urn() * (Lambda - 1.0 / Lambda) * k;
      double Ratio = knew / k;
      for (int i = 0; i < D; i++)
        NewMom[i] = OldMom[i] * Ratio;
      if (D == 2)
        Prop = 1.0;
      else if (D == 3)
        Prop = Ratio;
    }

    // if (isnan(Var.LoopMom[0][0]) || isnan(Var.LoopMom[0][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[1][0]) || isnan(Var.LoopMom[1][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[2][0]) || isnan(Var.LoopMom[2][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[3][0]) || isnan(Var.LoopMom[3][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[4][0]) || isnan(Var.LoopMom[4][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[5][0]) || isnan(Var.LoopMom[5][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[6][0]) || isnan(Var.LoopMom[6][0])) {
    //   cout << "Na" << endl;
    // }
  } else {
    for (int i = 0; i < D; i++)
      NewMom[i] = -OldMom[i];
    Prop = 1.0;
  }

  return Prop;
};

double markov::ShiftExtTransferK(const int &OldExtMomBin, int &NewExtMomBin) {
  NewExtMomBin = Random.irn(0, ExtMomBinSize - 1);
  return 1.0;
};

double markov::ShiftExtLegK(const momentum &OldExtMom, momentum &NewExtMom) {
  // double Theta = Random.urn() * 1.0 * PI;
  // NewExtMom[0] = Para.Kf * cos(Theta);
  // NewExtMom[1] = Para.Kf * sin(Theta);
  // return 1.0;

  int NewKBin = Random.irn(0, AngBinSize - 1);

  double AngCos = ver::Index2Angle(NewKBin, AngBinSize);
  double theta = acos(AngCos);
  NewExtMom[0] = Para.Kf * cos(theta);
  NewExtMom[1] = Para.Kf * sin(theta);

  // ASSERT_ALLWAYS(diag::Angle2Index(cos(theta), AngBinSize) == NewKBin,
  //                "Not matched, " << NewKBin << " vs "
  //                                << diag::Angle2Index(cos(theta),
  //                                AngBinSize));

  return 1.0;
};

double markov::ShiftTau(const double &OldTau, double &NewTau) {
  double x = Random.urn();
  if (x < 1.0 / 3) {
    double DeltaT = Para.Beta / 10.0;
    NewTau = OldTau + DeltaT * (Random.urn() - 0.5);
  } else if (x < 2.0 / 3) {
    NewTau = -OldTau;
  } else {
    NewTau = Random.urn() * Para.Beta;
  }

  // NewTau = Random.urn() * Para.Beta;
  if (NewTau < 0.0)
    NewTau += Para.Beta;
  if (NewTau > Para.Beta)
    NewTau -= Para.Beta;
  return 1.0;
}
