#include "weight.h"
#include "diagram.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

void weight::Initialization() {
  array<momentum *, 4> ExtLegK;

  if (DiagType == GAMMA) {
    ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1], &Var.LoopMom[2],
               &Var.LoopMom[3]};

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
      VerDiag.ResetMomMap(Ver4Root[order], ExtLegK);
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
    for (int order = 1; order <= Para.Order; order++) {
      Sigma[order] = BuildSigma(order + 1, &Var.LoopMom[0], &Var.LoopMom[1]);
      if (order < 4) {
        LOG_INFO(fmt::format(" Sigma, Order {0}\n", order));
        LOG_INFO(VerDiag.ToString(Sigma[order].Vertex));
      }
    }
  } else if (DiagType == POLAR) {

    /////////////////////////// Polar /////////////////////////
    ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1], &Var.LoopMom[2],
               &Var.LoopMom[3]};
    for (int order = 0; order <= Para.Order; order++) {
      Polar[order] = BuildPolar(order + 1, ExtLegK);
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
    Vertex4(Root, true);
  }
  double Factor = 1.0 / (Var.CurrAbsWeight * Para.ReWeight[Var.CurrOrder]);
  VerQTheta.Measure(Var.LoopMom[INL], Var.LoopMom[INR], Var.CurrExtMomBin,
                    Var.CurrOrder, ChanWeight, Factor);
}

void weight::Ver0(ver4 &Ver4) {
  // only bare coupling
  Ver4.Weight[0] = VerQTheta.Interaction(Ver4.LegK, 0.0, false, Ver4.InBox);
  // if (abs(Ver4.Weight[0][DIR] + 8.377) > 0.1) {
  //   cout << "Ver0 " << Ver4.Weight[0][DIR] << endl;
  //   ABORT("Err");
  // }
  // if (Ver4.Level == 1) {
  //   cout << "ver0: " << Ver4.Weight[0][DIR] << Ver4.Weight[0][EX] << endl;
  //   cout << "ver0 LegK0: " << (*Ver4.LegK[0]).norm() << endl;
  //   cout << "ver0 LegK1: " << (*Ver4.LegK[1]).norm() << endl;
  //   cout << "ver0 LegK2: " << (*Ver4.LegK[2]).norm() << endl;
  //   cout << "ver0 LegK3: " << (*Ver4.LegK[3]).norm() << endl;
  //   cout << Ver4.LegK[0] << ", " << Ver4.LegK[1] << ", " << Ver4.LegK[2] <<
  //   ", "
  //        << Ver4.LegK[3] << endl;
  // }
  return;
}

void weight::Benchmark(int LoopNum, diagram Diagram, int Step) {
  timer Timer;
  Timer.start();
  for (int i = 0; i < Step; i++)
    Evaluate(LoopNum);
  LOG_INFO(Timer << "s per " << Step << " step for LoopNum " << LoopNum
                 << " Diagram " << Diagram);
  return;
}