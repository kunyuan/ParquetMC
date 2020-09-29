#include "weight.h"
#include "diag/vertex4.h"
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
    // vector<channel> Chan = {I, T, U, S};
    // vector<channel> Chan = {T};
    // vector<channel> Chan = {U, UC};
    // vector<channel> Chan = {T};
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO("Generating order " << order);
      Gamma[order].Build(0,     // level
                         order, // loopNum
                         4,     // loop index of the first internal K
                         0,     // tau index of the InTL leg
                         Chan, RIGHT);
      if (order < 4)
        LOG_INFO(Gamma[order].ToString());
    }

  } else if (DiagType == SIGMA) {
    /////////////////////////// Sigma /////////////////////////
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO("Generating order " << order);
      Sigma[order].Build(order);
      if (order < 4)
        LOG_INFO(Sigma[order].ToString());
    }
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
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO("Generating order " << order);
      Delta[order].Build(order);
      if (order < 4)
        LOG_INFO(Delta[order].Vertex.ToString());
    }
  }
}

double weight::Evaluate(int Order) {
  if (Order == 0)
    return 1.0;

  // higher order
  if (DiagType == diagtype::GAMMA) {
    Gamma[Order].Evaluate(Var.LoopMom[INL],  // KInL
                          Var.LoopMom[OUTL], // KOutL
                          Var.LoopMom[INR],  // KInR
                          Var.LoopMom[OUTR], // KOutR
                          true);
    auto &ChanWeight = Gamma[Order].ChanWeight;
    for (auto &w : ChanWeight)
      // collapse all channel to I
      ChanWeight[0] += w;
    // cout << Order << ", " << ChanWeight[0][DIR] << ", " << ChanWeight[0][EX]
    //      << endl;
    return ChanWeight[0][DIR] + ChanWeight[0][EX] / SPIN;

  } else if (DiagType == diagtype::POLAR) {
    // polarization diagram
    double Weight = Polar[Order].Evaluate();
    return Weight;
  } else if (DiagType == diagtype::SIGMA) {
    // self-energy diagram
    double Weight = Sigma[Order].Evaluate();
    return Weight;
  } else if (DiagType == diagtype::DELTA) {
    // self-energy diagram
    double Weight = Delta[Order].Evaluate();
    return Weight;
  }
}

void weight::Measure() {
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);

  if (DiagType == diagtype::GAMMA) {
    if (Var.CurrOrder == 0) {
      // order zero
      GammaObs.Measure0(Factor);
    } else {
      Gamma[Var.CurrOrder].Evaluate(Var.LoopMom[INL], Var.LoopMom[OUTL],
                                    Var.LoopMom[INR], Var.LoopMom[OUTR], true);

      auto &ChanWeight = Gamma[Var.CurrOrder].ChanWeight;

      // double CosAng = diag::Angle3D(Var.LoopMom[INL], Var.LoopMom[INR]);
      // int AngleIndex = diag::Angle2Index(CosAng, Para.AngBinSize);
      // cout << AngleIndex << " vs " << Var.CurrExtAngBin << endl;

      GammaObs.Measure(Var.CurrOrder, Var.CurrExtMomBin, Var.CurrExtAngBin,
                       ChanWeight, Factor);
    }
  } else {
    Factor /= Para.TauGrid.weight[Var.CurrExtTauBin];
    // Polar, Sigma, Delta can be handled together
    if (Var.CurrOrder == 0)
      OneBodyObs.Measure0(Factor);
    else
      OneBodyObs.Measure(Var.CurrOrder, Var.CurrExtMomBin, Var.CurrExtTauBin,
                         Evaluate(Var.CurrOrder), Factor);
  }
}

void weight::SaveToFile() {
  if (DiagType == GAMMA)
    GammaObs.Save();
  else
    OneBodyObs.Save();
}

void weight::Test() {
  // cout << "start testing ..." << endl;
  if (DiagType == GAMMA)
    Gamma[Var.CurrOrder].Test();
  else if (DiagType == SIGMA)
    Sigma[Var.CurrOrder].Test();
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