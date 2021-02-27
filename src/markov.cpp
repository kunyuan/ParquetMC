//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#define FMT_HEADER_ONLY
#include "markov.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

extern parameter Para;
extern RandomFactory Random;
extern variable Var;

using namespace mc;
using namespace diag;
using namespace std;

int LastInterTauIdx(int Order) {
  // the last internal tau index
  if (DiagType == GAMMA)
    return Order;
  else if (DiagType == SIGMA)
    return Order - 2;
  else if (DiagType == POLAR)
    return Order - 1;
  else if (DiagType == DELTA)
    return Order - 2;
  else
    ASSERT(false, "Not Implemented!");
}

int FirstInterLoopIdx() {
  if (DiagType == GAMMA)
    return 4;
  else if (DiagType == SIGMA)
    return 1;
  else if (DiagType == POLAR)
    return 1;
  else if (DiagType == DELTA)
    return 1;
  else
    ASSERT(false, "Not Implemented!");
}

int LastInterLoopIdx(int Order) {
  // make sure Order>=1
  return FirstInterLoopIdx() + Order - 1;
}

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
    double NewTau;

    // Generate New Tau
    int NewTauIndex = LastInterTauIdx(Var.CurrOrder) + 1;
    if (NewTauIndex < 1) {
      Prop = 1.0;
    } else {
      Prop = GetNewTau(NewTau);
      Var.Tau[NewTauIndex] = NewTau;
    }

    // Generate New Mom
    static momentum NewMom;
    Prop *= GetNewK(NewMom);
    Var.LoopMom[LastInterLoopIdx(Var.CurrOrder) + 1] = NewMom;

  } else {
    // decrese order
    if (Var.CurrOrder == 0)
      return;

    Name = DECREASE_ORDER;
    NewOrder = Var.CurrOrder - 1;

    // Remove OldTau
    int TauToRemove = LastInterTauIdx(Var.CurrOrder);

    if (TauToRemove < 1)
      Prop = 1.0;
    else
      Prop = RemoveOldTau(Var.Tau[TauToRemove]);

    // Remove OldMom
    int LoopToRemove = LastInterLoopIdx(Var.CurrOrder);

    Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);
  }
  Proposed[Name][Var.CurrOrder] += 1;

  // Weight.ChangeGroup(NewGroup);
  NewAbsWeight = fabs(Weight.Evaluate(NewOrder));
  double R = Prop * NewAbsWeight * Para.ReWeight[NewOrder] / Var.CurrAbsWeight /
             Para.ReWeight[Var.CurrOrder];

  if (Random.urn() < R) {
    Accepted[Name][Var.CurrOrder]++;
    Var.CurrOrder = NewOrder;
    Var.CurrAbsWeight = NewAbsWeight;
  }
  return;
};

void markov::ChangeExtTau() {
  // Gamma diagram doesn't have ExtTau variable
  if (DiagType == GAMMA)
    return;

  Proposed[CHANGE_EXTTAU][Var.CurrOrder]++;

  int ExtTauIdx = MaxTauNum - 1;
  int OldTauBin = Var.CurrExtTauBin;
  double Prop = ShiftExtTau(OldTauBin, Var.CurrExtTauBin);

  Var.Tau[ExtTauIdx] = Para.TauGrid.grid[Var.CurrExtTauBin];

  NewAbsWeight = fabs(Weight.Evaluate(Var.CurrOrder));

  double R = Prop * NewAbsWeight * Para.TauGrid.weight[Var.CurrExtTauBin] /
             Var.CurrAbsWeight / Para.TauGrid.weight[OldTauBin];

  if (Random.urn() < R) {
    Accepted[CHANGE_EXTTAU][Var.CurrOrder]++;
    Var.CurrAbsWeight = NewAbsWeight;
  } else {
    Var.CurrExtTauBin = OldTauBin;
    Var.Tau[ExtTauIdx] = Para.TauGrid.grid[OldTauBin];
  }
};

void markov::ChangeTau() {

  if (LastInterTauIdx(Var.CurrOrder) < 1)
    return;

  int TauIndex = Random.irn(1, LastInterTauIdx(Var.CurrOrder));

  // int TauIndex = Random.irn(0, LastInterTauIdx(Var.CurrOrder));

  Proposed[CHANGE_TAU][Var.CurrOrder]++;

  double CurrTau = Var.Tau[TauIndex];
  double NewTau;
  double Prop = ShiftTau(CurrTau, NewTau);

  Var.Tau[TauIndex] = NewTau;

  NewAbsWeight = fabs(Weight.Evaluate(Var.CurrOrder));

  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_TAU][Var.CurrOrder]++;
    Var.CurrAbsWeight = NewAbsWeight;
  } else {
    Var.Tau[TauIndex] = CurrTau;
  }
};

void markov::ChangeMomentum() {
  if (Var.CurrOrder == 0)
    return;
  double Prop;
  static momentum CurrMom;

  int LoopIndex =
      Random.irn(FirstInterLoopIdx(), LastInterLoopIdx(Var.CurrOrder));

  // cout << Var.CurrOrder << ", " << FirstInterLoopIdx() << ", "
  //      << LastInterLoopIdx(Var.CurrOrder) << ", " << LoopIndex << endl;

  CurrMom = Var.LoopMom[LoopIndex];
  Prop = ShiftK(CurrMom, Var.LoopMom[LoopIndex]);

  Proposed[CHANGE_MOM][Var.CurrOrder]++;

  NewAbsWeight = fabs(Weight.Evaluate(Var.CurrOrder));
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_MOM][Var.CurrOrder]++;
    Var.CurrAbsWeight = NewAbsWeight;
  } else {
    Var.LoopMom[LoopIndex] = CurrMom;
  }
};

void markov::ChangeExtMomentum() {
  double Prop;
  int OldExtMomBin, OldAngBin;
  static momentum OldMom;

  if (DiagType == GAMMA) {
    // the INL, OUTL, OUTR momentum are fixed
    OldMom = Var.LoopMom[INR];
    OldAngBin = Var.CurrExtAngBin;
    Prop = ShiftExtLegK(OldAngBin, Var.CurrExtAngBin);

    double theta = acos(Para.AngleGrid.grid[Var.CurrExtAngBin]);
    Var.LoopMom[INR][0] = Para.Kf * cos(theta);
    Var.LoopMom[INR][1] = Para.Kf * sin(theta);
    Var.LoopMom[OUTR] = Var.LoopMom[INR];

  } else if (DiagType == SIGMA || DiagType == POLAR || DiagType == DELTA) {
    // In Momentum
    OldExtMomBin = Var.CurrExtMomBin;
    Prop = ShiftExtTransferK(OldExtMomBin, Var.CurrExtMomBin);
    if (DiagType == POLAR)
      Var.LoopMom[0][0] = Para.BoseKGrid.grid[Var.CurrExtMomBin];
    else
      Var.LoopMom[0][0] = Para.FermiKGrid.grid[Var.CurrExtMomBin];
  }

  Proposed[CHANGE_EXTMOM][Var.CurrOrder]++;

  NewAbsWeight = fabs(Weight.Evaluate(Var.CurrOrder));
  double R = Prop * NewAbsWeight / Var.CurrAbsWeight;

  if (Random.urn() < R) {
    Accepted[CHANGE_EXTMOM][Var.CurrOrder]++;
    Var.CurrAbsWeight = NewAbsWeight;
  } else {
    if (DiagType == GAMMA) {
      Var.CurrExtAngBin = OldAngBin;
      Var.LoopMom[INR] = OldMom;
      Var.LoopMom[OUTR] = OldMom;
    } else if (DiagType == SIGMA || DiagType == POLAR || DiagType == DELTA) {
      Var.CurrExtMomBin = OldExtMomBin;
      if (DiagType == POLAR)
        Var.LoopMom[0][0] = Para.BoseKGrid.grid[Var.CurrExtMomBin];
      else
        Var.LoopMom[0][0] = Para.FermiKGrid.grid[Var.CurrExtMomBin];
    }
  }
};

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
  double phi = 2.0 * PI * Random.urn();
  if (D == 3) {
    double theta = PI * Random.urn();
    if (theta == 0.0)
      return 0.0;
    NewMom[0] = KAmp * sin(theta) * cos(phi);
    NewMom[1] = KAmp * sin(theta) * sin(phi);
    NewMom[D - 1] = KAmp * cos(theta);
    return 2.0 * dK                    // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI                  // prop density of Phi
           * PI                        // prop density of Theta
           * sin(theta) * KAmp * KAmp; // Jacobian
  } else if (D == 2) {
    NewMom[0] = KAmp * cos(phi);
    NewMom[1] = KAmp * sin(phi);
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
    if (k < 1.0e-12) {
      Prop = 0.0;
      NewMom = OldMom;
    } else {
      const double Lambda = 1.5;
      double Ratio = 1.0 / Lambda + Random.urn() * (Lambda - 1.0 / Lambda);
      for (int i = 0; i < D; i++)
        NewMom[i] = OldMom[i] * Ratio;
      if (D == 2)
        Prop = 1.0;
      else if (D == 3)
        Prop = Ratio;
    }
  } else {
    for (int i = 0; i < D; i++)
      NewMom[i] = -OldMom[i];
    Prop = 1.0;
  }
  return Prop;
};

double markov::ShiftExtTransferK(const int &OldExtMomBin, int &NewExtMomBin) {
  if (DiagType == POLAR)
    NewExtMomBin = Random.irn(0, Para.BoseKGrid.size - 1);
  else
    NewExtMomBin = Random.irn(0, Para.FermiKGrid.size - 1);
  return 1.0;
};

double markov::ShiftExtTau(const int &OldTauBin, int &NewTauBin) {
  NewTauBin = Random.irn(0, Para.TauGrid.size - 1);
  return 1.0;
}

double markov::ShiftExtLegK(const int &OldExtMom, int &NewExtMom) {
  // double Theta = Random.urn() * 1.0 * PI;
  // NewExtMom[0] = Para.Kf * cos(Theta);
  // NewExtMom[1] = Para.Kf * sin(Theta);
  // return 1.0;

  NewExtMom = Random.irn(0, Para.AngleGrid.size - 1);

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

#define NAME(x) #x

markov::markov() {
  ///==== initialize Weight ============================//
  Weight.Initialization();

  //===== initialize updates related variable ==========//

  UpdatesName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdatesName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdatesName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdatesName[CHANGE_EXTMOM] = NAME(CHANGE_EXTMOM);
  UpdatesName[CHANGE_TAU] = NAME(CHANGE_TAU);
  UpdatesName[CHANGE_EXTTAU] = NAME(CHANGE_EXTTAU);

  // for(int i=0;i<MCUpdates;i++)
  // UpdatesName[(Updates)i]=NAME((Updates))

  InitialArray(&Accepted[0][0], 1.0e-10, MCUpdates * (MaxOrder + 1));
  InitialArray(&Proposed[0][0], 1.0e-10, MCUpdates * (MaxOrder + 1));

  AdjustGroupReWeight();
};

void markov::AdjustGroupReWeight(){};

std::string markov::_DetailBalanceStr(Updates op) {
  string Output = string(80, '-') + "\n";
  Output += UpdatesName[op] + ":\n";
  double TotalProposed = 0.0, TotalAccepted = 0.0;
  for (int i = 0; i <= Para.Order; i++) {
    if (!IsEqual(Proposed[op][i], 0.0)) {
      TotalAccepted += Accepted[op][i];
      TotalProposed += Proposed[op][i];
      Output +=
          fmt::sprintf("\t%8s%2i:%15g%15g%15g\n", "Order ", i, Proposed[op][i],
                       Accepted[op][i], Accepted[op][i] / Proposed[op][i]);
    }
  }
  if (!IsEqual(TotalProposed, 0.0)) {
    Output += fmt::sprintf("\t%10s:%15g%15g%15g\n", "Summation", TotalProposed,
                           TotalAccepted, TotalAccepted / TotalProposed);
  } else
    Output += "\tNo updates are proposed/accepted!\n";
  return Output;
}

void markov::PrintMCInfo() {
  string Output = "";
  Output = string(80, '=') + "\n";
  Output += "MC Counter: " + to_string(Var.Counter) + "\n";
  for (int i = 0; i < MCUpdates; i++)
    Output += _DetailBalanceStr((Updates)i);
  Output += string(80, '=') + "\n";
  LOG_INFO(Output);
}

void markov::PrintDeBugMCInfo() {
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Var.Counter) + ":\n";

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int d = 0; d < D; d++) {
    for (int i = 0; i < Var.CurrOrder + 3; i++)
      msg += ToString(Var.LoopMom[i][d]) + ", ";
    msg += "\n";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "Tau: \n";
  for (int i = 0; i < Var.CurrOrder + 1; i++)
    msg += ToString(Var.Tau[i]) + ", ";

  msg += "\n";

  LOG_INFO(msg);
}
