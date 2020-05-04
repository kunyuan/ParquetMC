#include "weight.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include <iostream>
#include <string>

using namespace diag;
using namespace std;

extern RandomFactory Random;
extern parameter Para;        // global parameters
extern variable Var;          // global MC variables
extern diag::propagator Prop; // global progator

void weight::Initialization() {
  array<momentum *, 4> ExtLegK;

  if (DiagType == GAMMA) {

    vector<channel> Chan = {I, T, U, S, TC, UC};
    // vector<channel> Chan = {T, TC};
    // vector<channel> Chan = {U, UC};
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO("Generating order " << order);
      Gamma[order].Build(0,     // level
                         order, // loopNum
                         4,     // loop index of the first internal K
                         0,     // tau index of the InTL leg
                         Chan, RIGHT, false);
      if (order < 4)
        LOG_INFO(Gamma[order].ToString());
    }

  } else if (DiagType == SIGMA) {
    /////////////////////////// Sigma /////////////////////////
  } else if (DiagType == POLAR) {
    /////////////////////////// Polar /////////////////////////
  } else if (DiagType == DELTA) {
    ////////////////////////// Delta /////////////////////////
  }
}

double weight::Evaluate(int LoopNum) {
  if (LoopNum == 0)
    return 1.0;

  // higher order
  if (DiagType == diagram::GAMMA) {
    Gamma[LoopNum].Evaluate(Var.LoopMom[INL],  // KInL
                            Var.LoopMom[OUTL], // KOutL
                            Var.LoopMom[INR],  // KInR
                            Var.LoopMom[OUTR], // KOutR
                            true);
    auto &ChanWeight = Gamma[LoopNum].ChanWeight;
    for (auto &w : ChanWeight)
      // collapse all channel to I
      ChanWeight[0] += w;
    return ChanWeight[0][DIR] + ChanWeight[0][EX] / SPIN;
  }
}

void weight::Measure() {
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);

  if (DiagType == GAMMA) {
    if (Var.CurrOrder != 0) {
      Gamma[Var.CurrOrder].Evaluate(Var.LoopMom[INL], Var.LoopMom[OUTL],
                                    Var.LoopMom[INR], Var.LoopMom[OUTR], true);

      auto &ChanWeight = Gamma[Var.CurrOrder].ChanWeight;

      GammaObs.Measure(Var.LoopMom[INL], Var.LoopMom[INR], Var.CurrExtMomBin,
                       Var.CurrOrder, ChanWeight, Factor);
    }
  }
}

void weight::SaveToFile() {
  if (DiagType == GAMMA)
    GammaObs.Save();
}

// void weight::MeasureSigma() {
//   double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
//   if (Var.CurrOrder == 0)
//     SigData.Measure0(Factor);
//   else if (Var.CurrOrder == 1) {
//     double Weight = EvaluateSigma(1, false);
//     SigData.Measure1(Var.CurrExtMomBin, Weight, Factor);
//   } else {
//     sigma &Root = Sigma[Var.CurrOrder];
//     EvaluateSigma(Var.CurrOrder, false);
//     // cout << Root.T[0] << "=" << Root.Weight[0] << ", " << Root.T[1] << "="
//     //      << Root.Weight[1] << endl;
//     SigData.Measure(Var.CurrOrder, Var.CurrExtMomBin, Root.T, Root.Weight,
//                     Factor);
//   }
// }
// void weight::MeasurePolar() {
//   double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
//   double Weight = EvaluatePolar(Var.CurrOrder);
//   PolarData.Measure(Var.CurrOrder, Var.CurrExtMomBin, Var.Tau[1] -
//   Var.Tau[0],
//                     Weight, Factor);
// }

// void weight::MeasureDelta() {
//   double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
//   double Weight = EvaluateDelta(Var.CurrOrder);
//   DeltaData.Measure(Var.CurrOrder, Var.CurrExtMomBin,
//                     Var.Tau[Var.CurrOrder - 1] - Var.Tau[0], Weight, Factor);
//   // Tau variable= the last Tau - the first tau
// }

void weight::Benchmark(int LoopNum, int Step) {
  timer Timer;
  Timer.start();
  for (int i = 0; i < Step; i++)
    Evaluate(LoopNum);
  LOG_INFO(Timer << "s per " << Step << " step for LoopNum " << LoopNum
                 << " Diagram " << DiagType);
  return;
}