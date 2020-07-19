#include "green.h"
#include <assert.h>
#include <math.h>

double fermiGreen(double beta, double tau, double Ek) {
  assert(tau > -beta && tau < beta);
  // equal time green's function
  if (tau == 0.0)
    tau = -1.0e-12;

  double green = 1.0;
  if (tau < 0.0) {
    tau += beta;
    green *= -1.0;
  }

  double x = beta * Ek / 2.0;
  double y = 2.0 * tau / beta - 1.0;
  if (x > 100.0)
    green *= exp(-x * (y + 1.0));
  else if (x < -100.0)
    green *= exp(x * (1.0 - y));
  else
    green *= exp(-x * y) / (2.0 * cosh(x));

  assert(isfinite(green));
  return green;
}

double fockYukawa(double k, double kF, double mass, bool shift) {
  // warning: this function only works for T=0!!!!
  double l = mass;
  double fock = 1.0 + l / kF * atan((k - kF) / l);
  fock -= l / kF * atan((k + kF) / l);
  fock -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
          log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
  fock *= (-2.0 * kF) / M_PI;

  if (shift) {
    double shift = 1.0 - l / kF * atan(2.0 * kF / l);
    shift -= l * l / 4.0 / kF / kF * log(l * l / (l * l + 4.0 * kF * kF));
    shift *= (-2.0 * kF) / M_PI;
    fock -= shift;
  }

  assert(isfinite(fock));
  return fock;
  // return fock;
}