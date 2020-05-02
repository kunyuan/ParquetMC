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
extern variable Var;

using namespace mc;
using namespace diag;
using namespace std;

#define NAME(x) #x
// #define COPYFROMTO(x, y)                                                       \
//   for (int i = 0; i < D; i++)                                                  \
//     y[i] = x[i];

// x=2*i, then PAIR(x)=2*i+1
// x=2*i+1, then PAIR(x)=2*i
#define PAIR(x) (int(x / 2) * 4 + 1 - x)

markov::markov() {
  ///==== initialize Weight ============================//
  Weight.Initialization();

  LOG_INFO("Initializating MC variables ...")
  // initialize momentum variables
  for (auto &mom : Var.LoopMom)
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn() * Para.Kf / sqrt(D);

  for (int i = 0; i < MaxTauNum; i++) {
    Var.Tau[i] = Random.urn() * Para.Beta;
  }

  // initialize spin variables
  // for (auto &sp : Var.LoopSpin)
  //   sp = (spin)(Random.irn(0, 1));

  // Var.LoopMom[0].fill(0.0);
  // Var.LoopMom[0] = Para.ExtMomTable[Var.CurrExtMomBin];

  if (DiagType == GAMMA) {
    Var.CurrExtMomBin = 0;
    for (int i = 1; i < D; i++) {
      Var.LoopMom[INL][i] = 0.0;
      Var.LoopMom[OUTL][i] = 0.0;
      Var.LoopMom[INR][i] = 0.0;
      Var.LoopMom[OUTR][i] = 0.0;
    }
    Var.LoopMom[INL][0] = Para.Kf;
    Var.LoopMom[OUTL][0] = Para.Kf;
    Var.LoopMom[INR][0] = Para.Kf;
    Var.LoopMom[OUTR][0] = Para.Kf;
  } else if (DiagType == SIGMA || DiagType == POLAR) {
    Var.CurrExtMomBin = ExtMomBinSize / 2;
    Var.LoopMom[0] = Para.ExtMomTable[Var.CurrExtMomBin];
  }

  Var.CurrTau = Var.Tau[1] - Var.Tau[0];

  // initialize group

  Var.Counter = 0;

  // Var.CurrGroup = &Groups[0];
  Var.CurrOrder = 0;

  // initialize RG staff
  // Var.CurrScale = ScaleBinSize - 1;

  LOG_INFO("Calculating the weights of all objects...")

  Var.CurrAbsWeight = fabs(Weight.Evaluate(Var.CurrOrder));

  LOG_INFO("Initializating variables done.")

  //===== initialize updates related variable ==========//

  UpdatesName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdatesName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdatesName[CHANGE_GROUP] = NAME(CHANGE_GROUP);
  UpdatesName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdatesName[CHANGE_EXTMOM] = NAME(CHANGE_EXTMOM);
  UpdatesName[CHANGE_TAU] = NAME(CHANGE_TAU);
  UpdatesName[CHANGE_SCALE] = NAME(CHANGE_SCALE);

  // for(int i=0;i<MCUpdates;i++)
  // UpdatesName[(Updates)i]=NAME((Updates))

  InitialArray(&Accepted[0][0], 1.0e-10, MCUpdates * (MaxOrder + 1));
  InitialArray(&Proposed[0][0], 1.0e-10, MCUpdates * (MaxOrder + 1));

  ///==== initialize observable =======================//
  ///=== Do all kinds of test  =======================//
  ///==== Set Reweighting factor =====================//
  AdjustGroupReWeight();
};

void markov::Measure() {
  if (DiagType == GAMMA)
    Weight.MeasureUST();
  else if (DiagType == SIGMA)
    Weight.MeasureSigma();
  else if (DiagType == POLAR)
    Weight.MeasurePolar();
};

void markov::AdjustGroupReWeight(){};

void markov::LoadFile() {
  Weight.SigData.LoadWeight();
  // if (DiagType == GAMMA){
  //   Weight.VerQTheta.LoadWeight();
  // }
  // else if (DiagType == SIGMA)
  //   Weight.SigData.LoadWeight();
  // else if (DiagType == POLAR)
  //   Weight.SigData.LoadWeight();
};

void markov::SaveToFile(bool Simple) {
  if (DiagType == GAMMA)
    Weight.VerQTheta.Save(Simple);
  else if (DiagType == SIGMA)
    Weight.SigData.Save();
  else if (DiagType == POLAR)
    Weight.PolarData.Save();
};

void markov::ClearStatis() { Weight.VerQTheta.ClearStatis(); }

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
      // fmt::format("\t%8s%4s:%15g%15g%15g\n", "Group", Groups[i].Name,
      //             Proposed[op][i], Accepted[op][i],
      //             Accepted[op][i] / Proposed[op][i]);
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

int markov::DynamicTest() {}

void markov::PrintDeBugMCInfo() {
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Var.Counter) + ":\n";

  // msg += string(80, '=') + "\n";
  // msg += "GWeight: \n";
  // for (int i = 0; i < Var.CurrGroup->GNum; i++)
  //   msg += ToString(Var.CurrGroup->Diag[0].G[i]->Weight) + "; ";
  // msg += "\n";

  // msg += "VerWeight: \n";
  // for (int i = 0; i < Var.CurrGroup->Ver4Num; i++) {
  //   msg += ToString(Var.CurrGroup->Diag[0].Ver4[i]->Weight[0]) + ", ";
  //   msg += ToString(Var.CurrGroup->Diag[0].Ver4[i]->Weight[1]) + "; ";
  // }
  // msg += "\n";

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
