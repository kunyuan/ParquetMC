#include "weight.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include <iostream>
#include <string>

using namespace diag;
using namespace std;

extern parameter Para;        // global parameters
extern variable Var;          // global MC variables
extern diag::propagator Prop; // global progator

void weight::Initialization() {
  array<momentum *, 4> ExtLegK;

  if (DiagType == GAMMA) {

    vector<channel> Chan = {I, T, U, S, TC, UC};
    // vector<channel> Chan = {T};
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
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO("Generating order " << order);
      Polar[order].Build(order);
      if (order < 4)
        LOG_INFO(Polar[order].Vertex.ToString());
    }
  } else if (DiagType == DELTA) {
    ////////////////////////// Delta /////////////////////////
  }
}

double weight::Evaluate(int Order) {
  if (Order == 0)
    return 1.0;

  // higher order
  if (DiagType == diagram::GAMMA) {
    Gamma[Order].Evaluate(Var.LoopMom[INL],  // KInL
                          Var.LoopMom[OUTL], // KOutL
                          Var.LoopMom[INR],  // KInR
                          Var.LoopMom[OUTR], // KOutR
                          true);
    auto &ChanWeight = Gamma[Order].ChanWeight;
    for (auto &w : ChanWeight)
      // collapse all channel to I
      ChanWeight[0] += w;
    return ChanWeight[0][DIR] + ChanWeight[0][EX] / SPIN;

  } else if (DiagType == diagram::POLAR) {
    // polarization diagram
    double Weight = Polar[Order].Evaluate();
    // cout << Order << ", " << Weight << endl;
    return Weight;
  }
}

void weight::Measure() {
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);

  if (DiagType == diagram::GAMMA) {
    // measure Gamma
    if (Var.CurrOrder == 0) {
      // order zero
      GammaObs.Measure0(Factor);
    } else {
      Gamma[Var.CurrOrder].Evaluate(Var.LoopMom[INL], Var.LoopMom[OUTL],
                                    Var.LoopMom[INR], Var.LoopMom[OUTR], true);

      auto &ChanWeight = Gamma[Var.CurrOrder].ChanWeight;

      GammaObs.Measure(Var.LoopMom[INL], Var.LoopMom[INR], Var.CurrExtMomBin,
                       Var.CurrOrder, ChanWeight, Factor);
    }

  } else if (DiagType == diagram::POLAR) {
    // measure polarization
    if (Var.CurrOrder == 0)
      PolarObs.Measure0(Factor);
    else
      PolarObs.Measure(Var.CurrOrder, Var.CurrExtMomBin,
                       Var.Tau[1] - Var.Tau[0], Polar[Var.CurrOrder].Evaluate(),
                       Factor);
  }
}

void weight::SaveToFile() {
  if (DiagType == GAMMA)
    GammaObs.Save();
  else if (DiagType == POLAR)
    PolarObs.Save();
}

void weight::Benchmark(int LoopNum, int Step) {
  timer Timer;
  Timer.start();
  for (int i = 0; i < Step; i++)
    Evaluate(LoopNum);
  LOG_INFO(Timer << "s per " << Step << " step for LoopNum " << LoopNum
                 << " Diagram " << DiagType);
  return;
}