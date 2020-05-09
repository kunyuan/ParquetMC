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

  if (DiagType == POLAR)
    // polarization
    Para.ReWeight = {0.05, 1.0, 2.0, 5.0, 0.2, 0.05, 1.0, 1.0, 1.0, 1.0};
  else if (DiagType == SIGMA)
    // Sigma
    Para.ReWeight = {1.0, 3.0, 30.0, 1.0, 0.2, 0.05, 1.0, 1.0, 1.0, 1.0};
  else if (DiagType == GAMMA)
    // Gamma
    Para.ReWeight = {10.0, 0.8, 0.4, 0.1, 0.4, 0.4, 1.0, 1.0, 1.0, 1.0};
  else if (DiagType == DELTA)
    // Delta
    Para.ReWeight = {1.0, 0.8, 10.0, 1.0, 0.4, 0.4, 1.0, 1.0, 1.0, 1.0};
  else
    ABORT("Not implemented!");

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

  LOG_INFO("Inverse Temperature: " << Para.Beta << "\n"
                                   << "r_s: " << Para.Rs << "\n"
                                   << "Fermi Mom: " << Para.Kf << "\n"
                                   << "Fermi Energy: " << Para.Ef << "\n");

  Para.PrinterTimer = 10;
  Para.SaveFileTimer = 10;
  Para.ReweightTimer = 30;
  Para.MessageTimer = 10;

  // Load external variable tables
  // try {
  LOG_INFO("Loading grids ...");

  ifstream File;
  string line;
  File.open("grid.data", ios::in);
  if (File.is_open()) {
    File >> Para.TauBinSize;
    LOG_INFO("TauTable: " << Para.TauBinSize);
    getline(File, line);
    Para.ExtTauTable.clear();
    double bin;
    for (int t = 0; t < Para.TauBinSize; ++t) {
      File >> bin;
      Para.ExtTauTable.push_back(bin);
    }

    File >> Para.ExtMomBinSize;
    LOG_INFO("MomTable: " << Para.ExtMomBinSize);
    getline(File, line);
    Para.ExtMomTable.clear();
    for (int k = 0; k < Para.ExtMomBinSize; ++k) {
      momentum mom;
      mom.setZero();
      File >> mom[0];
      Para.ExtMomTable.push_back(mom);
    }

    File >> Para.AngBinSize;
    getline(File, line);
    LOG_INFO("Angle Table: " << Para.AngBinSize);
    Para.AngleTable.clear();
    for (int ang = 0; ang < Para.ExtMomBinSize; ++ang) {
      double angle;
      File >> angle;
      Para.AngleTable.push_back(angle);
    }
    cout << "\n";

    // for (int i = 0; i < Para.AngBinSize; i++) {
    //   Para.AngleTable[i] = diag::Index2Angle(i, Para.AngBinSize);
    // }

    File.close();
  } else
    ABORT("Can not load grid file! \n");
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

  // reference tau, it should not be updated
  Var.Tau[0] = 0.0;
  // Var.Tau[0] = Para.Beta / 2.0;

  // Set the potential ExtTauBin
  Var.CurrExtTauBin = 0;
  Var.Tau[MaxTauNum - 1] = Para.ExtTauTable[Var.CurrExtTauBin];
  // cout << "Tau: " << Var.Tau[MaxTauNum - 1] << endl;

  if (DiagType == GAMMA) {
    Var.CurrExtMomBin = 0;
    for (int i = 0; i < 4; ++i)
      Var.LoopMom[i].setZero();
    Var.LoopMom[INL][0] = Para.Kf;
    Var.LoopMom[OUTL][0] = Para.Kf;

    Var.CurrExtAngBin = 0;
    double theta = acos(Para.AngleTable[Var.CurrExtAngBin]);
    Var.LoopMom[INR][0] = Para.Kf * cos(theta);
    Var.LoopMom[INR][1] = Para.Kf * sin(theta);

    Var.LoopMom[OUTR] = Var.LoopMom[INR];

  } else if (DiagType == SIGMA || DiagType == POLAR || DiagType == DELTA) {
    Var.CurrExtMomBin = 0;
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
  //   Markov.Weight.Benchmark(order, 10000);
  // }
  // exit(0);
  //////////////////////////////////////////////////////

  LOG_INFO("Start simulation ...")
  int Block = 0;
  while (Block < Para.TotalStep) {
    Block++;

    for (int i = 0; i < 1000000; i++) {
      Var.Counter++;

      double x = Random.urn();
      if (x < 1.0 / 5.0) {
        Markov.ChangeOrder();
      } else if (x < 2.0 / 5.0) {
        Markov.ChangeMomentum();
      } else if (x < 3.0 / 5.0) {
        Markov.ChangeExtMomentum();
      } else if (x < 4.0 / 5.0) {
        Markov.ChangeTau();
      } else if (x < 5.0 / 5.0) {
        Markov.ChangeExtTau();
      }

      // cout << Var.Tau[0] << endl;
      // cout << Var.CurrExtTauBin << ": " << Var.Tau[0] << "-> "
      //      << Var.Tau[MaxTauNum - 1] << endl;
      // Markov.Weight.Test();
      // cout << Var.LoopMom[0] << endl;
      // cout << Var.Tau[MaxTauNum - 1] << endl;
      // exit(0);
      if (Var.CurrOrder == -1) {

        momentum K1p = Var.LoopMom[4];
        momentum K2p = Var.LoopMom[5] + Var.LoopMom[INR] - K1p;

        cout << Var.CurrAbsWeight << " vs "
             << Prop.Green(Var.Tau[1] - Var.Tau[0], Var.LoopMom[4], UP, 0) *
                    Prop.Green(Var.Tau[0] - Var.Tau[2], Var.LoopMom[4], UP, 0) *
                    Prop.Green(Var.Tau[2] - Var.Tau[1], Var.LoopMom[5], UP, 0) *
                    Prop.Green(Var.Tau[1] - Var.Tau[2], K2p, UP, 0) /
                    pow(2.0 * PI, 6) * SPIN
             << endl;

        // cout <<"K1="<< Prop.Green(Var.Tau[1] - Var.Tau[0], Var.LoopMom[4],
        // UP, 0)<<endl; cout <<"K1p="<< Prop.Green(Var.Tau[0] - Var.Tau[1],
        // Var.LoopMom[4], UP, 0)<<endl;
        //             Prop.Green(Var.Tau[0] - Var.Tau[1], Var.LoopMom[4], UP,
        //             0) * Prop.Green(Var.Tau[2] - Var.Tau[1], Var.LoopMom[5],
        //             UP, 0) * Prop.Green(Var.Tau[1] - Var.Tau[2], K2p, UP, 0)
        //             / pow(2.0 * PI, 6) * 1.0
        //      << endl;

        // exit(0);
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
