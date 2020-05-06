#ifndef FeynCalc_global_h
#define FeynCalc_global_h

#include "utility/utility.h"
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <math.h>
#include <vector>

using namespace Eigen;
using namespace std;

enum type { GU, GW, RG, PARQUET, BARE, VARIATIONAL, RENORMALIZED };
enum diagram { SIGMA, POLAR, GAMMA, DELTA };

// turn off all assert
// #define NDEBUG

///////////  Global Constants ////////////////////
// D=2 or D=3
const int D = 3;
// spin index
const int SPIN = 2;
// type of diagram
const diagram DiagType = SIGMA;
// type of calculation
const type CalcType = VARIATIONAL;
// number of q bins of the external momentum
const int ExtMomBinSize = 32;
// number of bins for the angle between InL and InR legs
const int AngBinSize = 64;
// number of tau bin
const int TauBinSize = 256;
// Max diagram order
const int MaxOrder = 9;
const int MaxTauNum = MaxOrder + 1;
const int MaxMomNum = MaxOrder + 3;
// MaxMomNum = get_power<2, MaxOrder + 1>::value * 128;

// momentum vector
typedef Matrix<double, D, 1> momentum;

// vertex4 weight, has Direct and Exchange components
typedef Vector2d verWeight;

/////////// Global Parameter ////////////////////
struct parameter {
  // physical parameters
  double Rs, Ef, Kf;    // r_s, fermi energy, fermi momentum,
  double Nf, Mu, Beta;  // chemical potential, inverse temperature
  double Mass2, Lambda; // screening length^2, shift
  double Charge2;       // screening length^2
  double MaxExtMom;     // the maximum external momentum

  // MC inputs
  int Order;
  int TotalStep;                // total steps of the Monte Carlo
  int Seed;                     // rng seed
  int PID;                      // ID of the job
  int Sweep;                    // how many MC steps between two measuring
  std::vector<double> ReWeight; // reweight factor for each group

  // others
  int PrinterTimer;  // how many seconds between to printing to screen
  int SaveFileTimer; // how many secondes between saving to file
  int MessageTimer;  // how many secondes between two checking for message
  int ReweightTimer; // how many secondes between two reweighting

  // external variable tables
  // external bosonic Momentum (transfer momentum)
  momentum ExtMomTable[ExtMomBinSize];
  // external fermionic Momentum (LegK momentum)
  momentum ExtLegKTable[AngBinSize];
  double AngleTable[AngBinSize];
  double ExtTauTable[TauBinSize];
};

struct variable {
  long long int Counter; // counter to save the current MC step

  int CurrOrder;
  int CurrExtMomBin; // current bin of the external momentum
  int CurrTau;
  double CurrAbsWeight; // current abs weight

  // interval variables
  array<momentum, MaxMomNum> LoopMom; // all momentum loop variables
  array<double, MaxTauNum> Tau;       // all tau variables
};

//////////   Generic Global Constants  /////////////////
const double TM32 = 1.0 / (pow(2.0, 32));
const double EPS = 1.0e-9;
const int MAXINT = 2147483647;
const int MININT = -2147483647;
const double PI = 3.1415926535897932384626433832795;
const double MACHEPS = 2.22044604925031E-16; // Macheps + 1.0 > 1.0
const double MAXREAL = 1.0e30;
const double MINREAL = -1.0e30;

enum spin { DOWN, UP };
#define FLIPSPIN(x) spin(1 - x)
// Spin DOWN: 0,  Spin UP:1

#define FLIP(x) (1 - x)

const int IN = 0, OUT = 1;
const int LEFT = 0, RIGHT = 1;
const int INL = 0, OUTL = 1, INR = 2, OUTR = 3;
const int DIRECT = 0, EXCHANGE = 1;
const int DIR = 0, EX = 1;

//////////////////////////////////////////////////////

#define FMT_HEADER_ONLY

#endif