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

#define NAME(x) #x
// #define COPYFROMTO(x, y)                                                       \
//   for (int i = 0; i < D; i++)                                                  \
//     y[i] = x[i];

// x=2*i, then PAIR(x)=2*i+1
// x=2*i+1, then PAIR(x)=2*i
#define PAIR(x) (int(x / 2) * 4 + 1 - x)

markov::markov() : Var(Weight.Var) {
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

  Var.CurrExtMomBin = 0;

  // Var.LoopMom[0].fill(0.0);
  // Var.LoopMom[0] = Para.ExtMomTable[Var.CurrExtMomBin];

  for (int i = 1; i < D; i++) {
    Var.LoopMom[INL][i] = 0.0;
    Var.LoopMom[OUTL][i] = 0.0;
  }
  Var.LoopMom[INL][0] = Para.Kf;
  Var.LoopMom[OUTL][0] = Para.Kf;

  Var.CurrTau = Var.Tau[1] - Var.Tau[0];

  // initialize group

  Var.CurrVersion = 0;

  // Var.CurrGroup = &Groups[0];
  Var.CurrOrder = 0;

  Var.CurrIRScaleBin = ScaleBinSize / 1.5;

  Var.CurrChannel = dse::T;

  // initialize RG staff
  // Var.CurrScale = ScaleBinSize - 1;
  Var.CurrScale = Para.Kf;

  LOG_INFO("Calculating the weights of all objects...")

  Var.CurrWeight = Weight.Evaluate(Var.CurrOrder, Var.CurrChannel);
  Var.CurrAbsWeight = abs(Var.CurrWeight.Sum());

  LOG_INFO("Initializating variables done.")

  //===== initialize updates related variable ==========//

  UpdatesName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdatesName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdatesName[CHANGE_GROUP] = NAME(CHANGE_GROUP);
  UpdatesName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdatesName[CHANGE_TAU] = NAME(CHANGE_TAU);
  UpdatesName[CHANGE_SCALE] = NAME(CHANGE_SCALE);
  UpdatesName[VER2VER] = NAME(VER2VER);
  UpdatesName[SIGMA2VER] = NAME(SIGMA2VER);
  UpdatesName[VER2SIGMA] = NAME(VER2SIGMA);

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
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder] *
                         Para.ReWeightChan[Var.CurrChannel]);
  int TauIndex = GetTauNum(Var.CurrOrder, Var.CurrChannel) - 1;

  Weight.VerQTheta.Measure(Var.LoopMom[1], Var.LoopMom[2], Var.CurrExtMomBin,
                           Var.CurrOrder, Var.Tau[TauIndex] - Var.Tau[0],
                           Var.CurrChannel, Var.CurrWeight, Factor);
};

void markov::AdjustGroupReWeight(){};

void markov::LoadFile() { Weight.VerQTheta.LoadWeight(); };

void markov::SaveToFile(bool Simple) { Weight.VerQTheta.Save(Simple); };

void markov::ClearStatis() { Weight.VerQTheta.ClearStatis(); }

std::string markov::_DetailBalanceStr(Updates op) {
  string Output = string(80, '-') + "\n";
  Output += UpdatesName[op] + ":\n";
  double TotalProposed = 0.0, TotalAccepted = 0.0;
  for (int i = 0; i <= Para.Order; i++) {
    if (!Equal(Proposed[op][i], 0.0)) {
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
  if (!Equal(TotalProposed, 0.0)) {
    Output += fmt::sprintf("\t%10s:%15g%15g%15g\n", "Summation", TotalProposed,
                           TotalAccepted, TotalAccepted / TotalProposed);
  } else
    Output += "\tNo updates are proposed/accepted!\n";
  return Output;
}

void markov::PrintMCInfo() {
  string Output = "";
  Output = string(80, '=') + "\n";
  Output += "MC Counter: " + to_string(Para.Counter) + "\n";
  for (int i = 0; i < MCUpdates; i++)
    Output += _DetailBalanceStr((Updates)i);
  Output += string(80, '=') + "\n";
  LOG_INFO(Output);
}

int markov::DynamicTest() {}

void markov::PrintDeBugMCInfo() {
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Para.Counter) + ":\n";

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
