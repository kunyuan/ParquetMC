#include "weight.h"
#include "diagram.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

extern variable Var;

void weight::Initialization() {
  array<momentum *, 4> ExtLegK;

  if (DiagType == GAMMA) {
    // ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1], &Var.LoopMom[2],
    //            &Var.LoopMom[3]};

    vector<channel> Chan = {I, T, U, S, TC, UC};
    // vector<channel> Chan = {T, TC};
    // vector<channel> Chan = {U, UC};
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO("Generating order " << order);
      Ver4Root[order] = VerDiag.Vertex(0,     // level
                                       order, // loopNum
                                       4, // loop index of the first internal K
                                       0, // tau index of the InTL leg
                                       Chan, RIGHT, false);
      // VerDiag.ResetMomMap(Ver4Root[order], ExtLegK);
      // the vertex LegK must be initialized after the memory allocaiton of the
      // vertex tree
      if (order < 4)
        LOG_INFO(VerDiag.ToString(Ver4Root[order]));
    }
    for (int order = 1; order <= Para.Order; order++) {
      LOG_INFO(fmt::format(" Order {0}: ExtT Num {1}\n", order,
                           Ver4Root[order].T.size()));
    }

  } else if (DiagType == SIGMA) {
    /////////////////////////// Sigma /////////////////////////
    for (int order = 2; order <= Para.Order; order++) {
      Sigma[order] = BuildSigma(order);
      // SetSigmaMom(Sigma[order], &Var.LoopMom[0], &Var.LoopMom[1]);

      if (order < 4) {
        LOG_INFO(fmt::format(" Sigma, Order {0}\n", order));
        LOG_INFO(VerDiag.ToString(Sigma[order].Vertex));
      }
    }
  } else if (DiagType == POLAR) {
    /////////////////////////// Polar /////////////////////////
    for (int order = 2; order <= Para.Order; order++) {
      Polar[order] = BuildPolar(order);
      // if (order < 4) {
      //   LOG_INFO(fmt::format(" Polar, Order {0}\n", order));
      //   LOG_INFO(VerDiag.ToString(Sigma[order].Vertex));
      // }
    }
  } else if (DiagType == DELTA) {
    for (int order = 1; order <= Para.Order; order++) {
      Delta[order] = BuildDelta(order);
      // if (order < 4) {
      //   LOG_INFO(fmt::format(" Polar, Order {0}\n", order));
      //   LOG_INFO(VerDiag.ToString(Sigma[order].Vertex));
      // }
    }
  }
}

void weight::MeasureUST() {
  if (Var.CurrOrder != 0) {
    ver4 &Root = Ver4Root[Var.CurrOrder];
    // if (Root.Weight.size() != 0) {
    Vertex4(Root, Var.LoopMom[INL], Var.LoopMom[OUTL], Var.LoopMom[INR],
            Var.LoopMom[OUTR], true);
  }
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
  VerQTheta.Measure(Var.LoopMom[INL], Var.LoopMom[INR], Var.CurrExtMomBin,
                    Var.CurrOrder, ChanWeight, Factor);
}

void weight::MeasureSigma() {
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
  if (Var.CurrOrder == 0)
    SigData.Measure0(Factor);
  else if (Var.CurrOrder == 1) {
    double Weight = EvaluateSigma(1, false);
    SigData.Measure1(Var.CurrExtMomBin, Weight, Factor);
  } else {
    sigma &Root = Sigma[Var.CurrOrder];
    EvaluateSigma(Var.CurrOrder, false);
    // cout << Root.T[0] << "=" << Root.Weight[0] << ", " << Root.T[1] << "="
    //      << Root.Weight[1] << endl;
    SigData.Measure(Var.CurrOrder, Var.CurrExtMomBin, Root.T, Root.Weight,
                    Factor);
  }
}
void weight::MeasurePolar() {
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
  double Weight = EvaluatePolar(Var.CurrOrder);
  PolarData.Measure(Var.CurrOrder, Var.CurrExtMomBin, Var.Tau[1] - Var.Tau[0],
                    Weight, Factor);
}

void weight::MeasureDelta() {
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
  double Weight = EvaluateDelta(Var.CurrOrder);
  DeltaData.Measure(Var.CurrOrder, Var.CurrExtMomBin,
                    Var.Tau[Var.CurrOrder - 1] - Var.Tau[0], Weight, Factor);
  // Tau variable= the last Tau - the first tau
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