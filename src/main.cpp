#include "grid.h"
#include "markov.h"
#include "utility/timer.h"
#include <algorithm>
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
void InitPara(), InitVar();

const string HelpStr = "Two parameters: PID Seed";

int main(int argc, const char *argv[]) {
#ifdef NDEBUG
  LOG_INFO("NDEBUG mode is OFF.");
#else
  LOG_INFO("NDEBUG mode is ON.");
#endif
  // take two parameters: PID and Seed
  Para.PID = atoi(argv[1]);
  Para.Seed = atoi(argv[2]);

  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

  ASSERT_ALLWAYS(Para.Seed > 0, "Random number seed must be positive integer!");
  ASSERT_ALLWAYS(Para.PID >= 0, "PID must be positive integer!");
  Random.Reset(Para.Seed);

  InitPara(); // initialize global parameters

  markov Markov;
  InterruptHandler Interrupt;

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

      // cout << Var.LoopMom[1].norm() << endl;

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

stringstream GetLine(ifstream &File) {
  string line;
  while (true) {
    getline(File, line);
    line = trim(line);
    // cout << "get " << line << endl;
    if (line.size() > 0 && line[0] != '#') {
      replace(line.begin(), line.end(), ',', ' ');
      // cout << "return " << line << endl;
      return stringstream(line);
    }
  }
}

void InitPara() {

  ifstream File;
  string line;
  File.open("parameter", ios::in);
  ASSERT_ALLWAYS(File.is_open(), "Can not load parameters! \n");
  // parameters
  auto paraStream = GetLine(File);
  paraStream >> Para.Order >> Para.Beta >> Para.Rs >> Para.Mass2 >>
      Para.Lambda >> Para.Charge2 >> Para.TotalStep;

  // grid information
  int RealFreqGridSize;
  double MaxRealFreq;
  auto gridStream = GetLine(File);
  gridStream >> Para.TauBinSize >> Para.ExtMomBinSize >> Para.AngBinSize >>
      RealFreqGridSize >> MaxRealFreq >> Para.TauBasisSize;

  // Timer information
  auto timerStream = GetLine(File);
  timerStream >> Para.PrinterTimer >> Para.SaveFileTimer >>
      Para.ReweightTimer >> Para.MessageTimer;

  // ReWeight information
  auto reweightStream = GetLine(File);
  for (int o = 0; o < Para.Order + 1; ++o)
    reweightStream >> Para.ReWeight[o];

  File.close();

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

  LOG_INFO("PrintTimer: " << Para.PrinterTimer << "\n"
                          << "SaveTimer: " << Para.SaveFileTimer << "\n"
                          << "ReWeightTimer: " << Para.ReweightTimer << "\n"
                          << "MessageTimer: " << Para.MessageTimer << "\n");

  // Load external variable tables
  // try {
  // LOG_INFO("Loading grids ...");

  // auto Grid = grid();
  // Grid.Initialize({0.0, Para.Beta / 2.0}, Para.TauBinSize / 2 + 1, true,
  //                 Para.Ef * Para.Beta / 4.0);

  int Size;
  File.open("grid.data", ios::in);
  ASSERT_ALLWAYS(File.is_open(), "Can not load grid file! \n");

  GetLine(File) >> Size;
  ASSERT_ALLWAYS(Size == Para.TauBinSize, "TauBinSize is invalid!");
  Para.ExtTauTable.clear();
  double bin;
  for (int t = 0; t < Para.TauBinSize; ++t) {
    File >> bin;
    Para.ExtTauTable.push_back(bin);
  }

  // auto TauGrid = tauGrid();
  // TauGrid.Initialize(Para.Beta, Para.TauBinSize, Para.Ef * Para.Beta / 10.0);
  // cout << TauGrid.ToString() << endl;
  // for (int t = 0; t < Para.TauBinSize; ++t) {
  //   Para.ExtTauTable[t] = TauGrid.Grid(t);
  // }

  GetLine(File) >> Size;
  ASSERT_ALLWAYS(Size == Para.TauBinSize, "TauBinSize is invalid!");
  Para.ExtTauReWeight.clear();
  for (int t = 0; t < Para.TauBinSize; ++t) {
    File >> bin;
    Para.ExtTauReWeight.push_back(bin);
  }

  // make the first and the last tau weights bigger
  // Para.ExtTauReWeight[0] /= 100.0;
  // Para.ExtTauReWeight[Para.TauBinSize - 1] /= 100.0;

  // normalize the tau bin weights
  double TotalWeight = 0.0;
  for (auto &w : Para.ExtTauReWeight)
    TotalWeight += w;
  for (auto &w : Para.ExtTauReWeight)
    w /= TotalWeight;

  GetLine(File) >> Size;
  ASSERT_ALLWAYS(Size == Para.ExtMomBinSize, "ExtMomBinSize is invalid!");
  Para.ExtMomTable.clear();
  for (int k = 0; k < Para.ExtMomBinSize; ++k) {
    momentum mom;
    mom.setZero();
    File >> mom[0];
    Para.ExtMomTable.push_back(mom);
  }

  GetLine(File) >> Size;
  ASSERT_ALLWAYS(Size == Para.AngBinSize, "AngBinSize is invalid!");
  Para.AngleTable.clear();
  for (int ang = 0; ang < Para.ExtMomBinSize; ++ang) {
    double angle;
    File >> angle;
    Para.AngleTable.push_back(angle);
  }
  File.close();
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
