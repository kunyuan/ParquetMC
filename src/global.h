#ifndef FeynCalc_global_h
#define FeynCalc_global_h

#include "utility/utility.h"
#include <Eigen/Dense>
#include <array>
#include <math.h>
#include <vector>

// turn off all assert
// #define NDEBUG

enum type { GU, GW, RG, PARQUET, BARE, VARIATIONAL, RENORMALIZED };
enum diagram { SIGMA, POLAR, GAMMA, DELTA };

///////////  Global Constants ////////////////////
const int D = 3;                   // dimensions, 2 or 3
const int SPIN = 2;                // spin index
const type CalcType = VARIATIONAL; // calculation type
const diagram DiagType = SIGMA;    // diagram type
const int ExtMomBinSize = 32;      // external K bins
const int AngBinSize = 64;         // angle bins
const int TauBinSize = 128;        // tau bin
const int MaxOrder = 9;            // Max diagram order
const int MaxTauNum = MaxOrder + 1;
const int MaxMomNum = MaxOrder + 3;
// MaxMomNum = get_power<2, MaxOrder + 1>::value * 128;

// momentum vector
typedef Eigen::Matrix<double, D, 1> momentum;
// vertex4 weight, has Direct and Exchange components
typedef Eigen::Vector2d verWeight;

/////////// Global Parameter ////////////////////
struct parameter {
  // physical parameters
  int Order;
  double Rs, Ef, Kf;    // r_s, fermi energy, fermi momentum,
  double Nf, Mu, Beta;  // chemical potential, inverse temperature
  double Mass2, Lambda; // screening length^2, shift
  double Charge2;       // screening length^2
  double MaxExtMom;     // the maximum external momentum

  // MC inputs
  int TotalStep;                // total steps of the Monte Carlo
  int Sweep;                    // how many MC steps between two measuring
  int Seed, PID;                // rng seed, job ID
  std::vector<double> ReWeight; // reweight factor for each group

  // others
  int PrinterTimer;  // time interval to print to screen
  int SaveFileTimer; // time interval to file
  int MessageTimer;  // time interval to check for new input data
  int ReweightTimer; // time interval to reweight different orders

  // external variable tables
  momentum ExtMomTable[ExtMomBinSize];
  momentum ExtLegKTable[AngBinSize];
  double AngleTable[AngBinSize];
  double ExtTauTable[TauBinSize];
};

struct variable {
  long long int Counter; // counter to save the current MC step

  int CurrOrder;        // current order
  double CurrAbsWeight; // current abs weight

  // external variables
  int CurrExtMomBin;
  int CurrExtTauBin;
  int CurrExtAngBin;

  // interval variables
  std::array<momentum, MaxMomNum> LoopMom; // all momentum loop variables
  std::array<double, MaxTauNum> Tau;       // all tau variables
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