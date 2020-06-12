#ifndef response_H
#define response_H

// return the bare fermionic Green's function
// Ek=k^2-mu
/**
 * @brief the bare Fermionic Green's function
 * @param Ek k^2/2m-μ
 * @param tau (-beta, beta)
 */
double fermiGreen(double beta, double tau, double Ek);

/**
 * @brief zero-temperature Fock diagram with a Yukawa type interaction
 *
 * @param k momentum
 * @param kF fermi momentum
 * @param mass Mass of the bosonic mode
 * @return fock diagram without a chemical-potential shift
 */

double fockYukawa(double k, double kF, double mass, bool shift);

#endif