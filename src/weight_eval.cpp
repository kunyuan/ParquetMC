#include "diagram.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include "utility/vector.h"
#include "weight.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

double weight::Evaluate(int LoopNum) {
  if (DiagType == diagram::GAMMA)
    return EvaluateGamma(LoopNum);
  else if (DiagType == diagram::SIGMA)
    return EvaluateSigma(LoopNum);
  else if (DiagType == diagram::POLAR)
    return EvaluatePolar(LoopNum);
  else if (DiagType == diagram::DELTA)
    return EvaluateDelta(LoopNum);
}

double weight::EvaluateGamma(int LoopNum) {
  if (LoopNum == 0) {
    // normalization
    return 1.0;
  } else {
    // if (Channel != dse::T)
    //   return 0.0;
    ver4 &Root = Ver4Root[LoopNum];
    if (Root.Weight.size() != 0) {
      Vertex4(Root,
              Var.LoopMom[INL],  // KInL
              Var.LoopMom[OUTL], // KOutL
              Var.LoopMom[INR],  // KInR
              Var.LoopMom[OUTR], // KOutR
              true);
      for (auto &w : ChanWeight)
        // collapse all channel to I
        ChanWeight[0] += w;
      // cout << ChanWeight[0].Abs() << endl;
      return ChanWeight[0].Abs();
    } else
      return 0.0;
  }
}

double weight::EvaluatePolar(int LoopNum) {
  double Factor = 1.0 / pow(2.0 * PI, D);
  // normalization
  if (LoopNum == 0)
    return 1.0;
  else if (LoopNum == 1) {
    double Tau = Var.Tau[1] - Var.Tau[0];
    double Weight = Fermi.Green(Tau, Var.LoopMom[1], UP, 0, 0);
    Weight *= Fermi.Green(-Tau, Var.LoopMom[1] - Var.LoopMom[0], UP, 0, 0);

    return -SPIN * Weight * Factor;
  }

  // loop order >=2
  dse::polar &P = Polar[LoopNum];
  ver4 &Ver4 = P.Vertex;
  momentum KOutL = Var.LoopMom[1] - Var.LoopMom[0];
  momentum KOutR = Var.LoopMom[2] + Var.LoopMom[0];

  Vertex4(Ver4,
          Var.LoopMom[1], // KInL
          KOutL,          // KOutL
          Var.LoopMom[2], // KInR
          KOutR,          // KOutR
          false);

  // cout << Ver4.Weight[0][DIR] << ", " << Ver4.Weight[0][EX] << endl;
  // cout << Var.LoopMom[0][0] << ", Mom1=" << Var.LoopMom[1].norm() << endl;
  // cout << "DirQ: " << (Var.LoopMom[1] - KOutL).norm() << endl;

  // evaluate all possible G
  EvaluateG(P.G[INL], Var.LoopMom[1]);
  EvaluateG(P.G[OUTL], KOutL);
  EvaluateG(P.G[INR], Var.LoopMom[2]);
  EvaluateG(P.G[OUTR], KOutR);

  int Size = Ver4.T.size();
  P.Weight = 0.0;
  for (int i = 0; i < Size; ++i) {
    auto &Gidx = P.Gidx[i];
    double Weight =
        (SPIN * SPIN * Ver4.Weight[i][DIR] + SPIN * Ver4.Weight[i][EX]);

    // attach four G
    for (int j = 0; j < 4; ++j) {
      Weight *= P.G[j][Gidx[j]].Weight;
      // cout << j << "=" << P.G[j][Gidx[j]].Weight << endl;
    }

    P.Weight += Weight;
  }
  // cout << P.Weight << endl;

  return P.Weight * Factor * Factor;
}

double weight::EvaluateSigma(int LoopNum, bool IsFast) {
  // normalization
  double Factor = 1.0 / pow(2.0 * PI, D);
  if (LoopNum == 0) {
    return 1.0;
  } else if (LoopNum == 1) {
    double GWeight =
        Fermi.Green(-EPS, Var.LoopMom[0] + Var.LoopMom[1], UP, 0, 0);
    double VerWeight = VerQTheta.Interaction(Var.LoopMom[1]);
    // cout << GWeight << ", " << VerWeight << endl;
    return GWeight * VerWeight * Factor;
  } else {
    // Sigma with LoopNum>=2
    dse::sigma &Si = Sigma[LoopNum];
    ver4 &Ver4 = Si.Vertex;

    Vertex4(Ver4, Var.LoopMom[0], // KInL
            Var.LoopMom[1],       // KOutL
            Var.LoopMom[1],       // KInR
            Var.LoopMom[0],       // KOutR
            false);

    EvaluateG(Si.G, Var.LoopMom[1]);
    if (IsFast)
      Si.Weight[0] = 0.0;
    else {
      for (auto &w : Si.Weight)
        w = 0.0;
    }
    int Size = Ver4.Weight.size();
    // Size = 1;
    for (int i = 0; i < Size; ++i) {

      auto &Weight = Ver4.Weight[i];
      int Gidx = Si.Gidx[i];
      double w = (Weight[DIR] + Weight[EX] * SPIN) * Si.G[Gidx].Weight;

      if (Si.T[i] != 0)
        w *= 0.5;

      if (IsFast)
        Si.Weight[0] += w * Factor;
      else
        Si.Weight[Si.SigTidx[i]] += w * Factor;
    }
    return Si.Weight[0];
  }
}

double weight::EvaluateDelta(int LoopNum) { return 0.0; }