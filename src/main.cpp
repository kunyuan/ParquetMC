#include "markov.h"
#include "utility/timer.h"
#include <iostream>
#include <math.h>

using namespace std;
using namespace mc;
void InitPara();
void InitVar();
void MonteCarlo();

// Global variable
RandomFactory Random;
parameter Para;        // global parameters
diag::propagator Prop; // global progator
variable Var;

void InitPara() {
  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

#ifdef NDEBUG
  LOG_INFO("NDEBUG mode is OFF.");
#else
  LOG_INFO("NDEBUG mode is ON.");
#endif

  // Para.ReWeight = {2.0, 0.8, 0.4, 0.4, 0.4, 0.4, 1.0, 1.0, 1.0, 1.0};
  Para.ReWeight = {1.0, 1.0, 1.0, 0.2, 0.2, 0.05, 1.0, 1.0, 1.0, 1.0};

  //// initialize the global parameter //////////////////////
  double Kf;
  if (D == 3) {
    Kf = pow(9.0 * PI / 4.0, 1.0 / 3.0) / Para.Rs; // 3D
  } else if (D == 2) {
    Kf = sqrt(2.0) / Para.Rs; // 2D
  } else {
    ABORT("Dimension " << D << " has not yet been implemented!");
  }
  Para.Kf = Kf;
  Para.Ef = Kf * Kf;
  Para.Mu = Para.Ef;
  Para.Nf = Kf / (4.0 * PI * PI) * SPIN;
  Para.MaxExtMom *= Kf;

  // scale all energy with E_F
  Para.Beta /= Para.Ef;

  for (int i = 0; i < AngBinSize; i++) {
    Para.AngleTable[i] = diag::Index2Angle(i, AngBinSize);
  }

  // initialize external momentum
  for (int i = 0; i < ExtMomBinSize; i++) {
    // the external momentum only has x component
    Para.ExtMomTable[i].setZero();
    if (DiagType == GAMMA)
      Para.ExtMomTable[i][0] = i * Para.MaxExtMom / ExtMomBinSize;
    else
      Para.ExtMomTable[i][0] = (i + 0.5) * Para.MaxExtMom / ExtMomBinSize;
  }

  for (int i = 0; i < TauBinSize; i++) {
    Para.ExtTauTable[i] = Para.Beta * i / TauBinSize + 1.0e-9;
  }
  // Para.ExtMomTable[0][0] = 0.0;
  // Para.ExtMomTable[1][0] = 2. * Para.Kf;

  LOG_INFO("Inverse Temperature: " << Para.Beta << "\n"
                                   << "r_s: " << Para.Rs << "\n"
                                   << "Fermi Mom: " << Para.Kf << "\n"
                                   << "Fermi Energy: " << Para.Ef << "\n");

  Para.PrinterTimer = 10;
  Para.SaveFileTimer = 10;
  Para.ReweightTimer = 30;
  Para.MessageTimer = 10;
}

void InitVar() {
  // initialize group
  Var.Counter = 0;
  // Var.CurrGroup = &Groups[0];
  Var.CurrOrder = 0;

  // initialize momentum variables
  for (auto &mom : Var.LoopMom)
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn() * Para.Kf / sqrt(D);

  for (auto &t : Var.Tau)
    t = Random.urn() * Para.Beta;
  Var.Tau[0] = 0.0; // reference tau

  if (DiagType == GAMMA) {
    Var.CurrExtMomBin = 0;
    for (int i = 0; i < 4; ++i) {
      Var.LoopMom[i].setZero();
      Var.LoopMom[i][0] = Para.Kf;
    }
  } else if (DiagType == SIGMA || DiagType == POLAR) {
    Var.CurrExtMomBin = ExtMomBinSize / 2;
    Var.LoopMom[0] = Para.ExtMomTable[Var.CurrExtMomBin];
  }
}

int main(int argc, const char *argv[]) {
  cout << "Order, Beta, Rs, Mass2, Lambda, Charge2, MaxExtMom(*kF), "
          "TotalStep(*1e6), "
          "Seed, "
          "PID\n";
  cin >> Para.Order >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.Lambda >>
      Para.Charge2 >> Para.MaxExtMom >> Para.TotalStep >> Para.Seed >> Para.PID;

  InitPara(); // initialize global parameters

  markov Markov;
  InterruptHandler Interrupt;

  Random.Reset(Para.Seed);

  timer ReweightTimer, PrinterTimer, SaveFileTimer, MessageTimer;
  PrinterTimer.start();
  SaveFileTimer.start();
  MessageTimer.start();
  ReweightTimer.start();

  LOG_INFO("Loading Weight ...")
  Markov.Weight.LoadFile();
  InitVar(); // initialize MC variables
  Var.CurrAbsWeight = fabs(Markov.Weight.Evaluate(Var.CurrOrder));

  ///////////////  Benchmark ////////////////////////////
  // for (int order = 1; order <= Para.Order; ++order) {
  //   Markov.Weight.Benchmark(order, 1000);
  // }
  //////////////////////////////////////////////////////

  LOG_INFO("Start simulation ...")
  int Block = 0;
  while (Block < Para.TotalStep) {
    Block++;

    for (int i = 0; i < 1000000; i++) {
      Var.Counter++;

      double x = Random.urn();
      if (x < 1.0 / 4.0) {
        Markov.ChangeOrder();
        // ;
      } else if (x < 2.0 / 4.0) {
        Markov.ChangeMomentum();
        // ;
      } else if (x < 3.0 / 4.0) {
        Markov.ChangeExtMomentum();
        // ;
      } else if (x < 4.0 / 4.0) {
        Markov.ChangeTau();
        // else if (x < 5.0 / 5.0) {
        //   Markov.ChangeScale();
        // ;
      }

      if (i % 8 == 0)
        // fast operations
        Markov.Weight.Measure();

      if (i % 1000 == 0) {
        // slow operations
        if (PrinterTimer.check(Para.PrinterTimer)) {
          Markov.Weight.Test();
          Markov.PrintDeBugMCInfo();
          Markov.PrintMCInfo();
          LOG_INFO(ProgressBar((double)Block / Para.TotalStep));
        }

        if (SaveFileTimer.check(Para.SaveFileTimer)) {
          Interrupt.Delay(); // the process can not be killed in saving
          Markov.Weight.SaveToFile();
          Interrupt.Resume(); // after this point, the process can be killed
        }

        if (ReweightTimer.check(Para.ReweightTimer)) {
          Markov.AdjustGroupReWeight();
          Para.ReweightTimer *= 1.5;
        }

        if (MessageTimer.check(Para.MessageTimer)) {
          LOG_INFO("Loading Weight...")
          Markov.Weight.LoadFile();
        }
      }
    }
  }

  Markov.PrintMCInfo();
  Interrupt.Delay(); // the process can not be killed in saving
  Markov.Weight.SaveToFile();
  Interrupt.Resume(); // after this point, the process can be killed

  LOG_INFO("Simulation is ended!");

  return 0;
}
