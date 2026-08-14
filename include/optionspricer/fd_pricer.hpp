#pragma once

#include "optionspricer/black_scholes.hpp"
#include "optionspricer/types.hpp"

namespace optionspricer {

// Crank-Nicolson finite-difference pricer for European and American options.
//
//   Black-Scholes PDE in backward time (tau = T - t):
//       V_tau = 0.5*sigma^2*S^2*V_SS + (r-q)*S*V_S - r*V
//
//   Uniform grid in S: S_i = i*dS, i = 0..M. Central differences give
//       (L V)_i = a_i*V_{i-1} + b_i*V_i + c_i*V_{i+1}
//   with dS cancelling entirely:
//       a_i = 0.5*sigma^2*i^2 - 0.5*(r-q)*i
//       b_i = -sigma^2*i^2 - r
//       c_i = 0.5*sigma^2*i^2 + 0.5*(r-q)*i
//
//   Theta-scheme:  (I - th*dt*L) V^{n+1} = (I + (1-th)*dt*L) V^n
//       th = 0.5 -> Crank-Nicolson (2nd order in time)
//       th = 1.0 -> fully implicit  (1st order, but strongly damping)
//
//   Rannacher smoothing: the first rannacher_steps time steps run fully
//   implicit to damp the oscillations Crank-Nicolson produces against the
//   kink in the payoff. Without it, gamma near the strike rings badly.
//
// Reuses BSResult for its output. price, delta, and gamma come off the
// grid (delta/gamma via central differences around the node S0 is pinned
// to); vega, theta, and rho are not computed by this pricer and are left
// at BSResult's default of 0.0.
//
//   S     - spot price
//   K     - strike price
//   r     - continuously compounded risk-free rate
//   q     - continuously compounded dividend yield
//   sigma - volatility (annualized)
//   T     - time to maturity in years
//   M     - number of spatial (asset price) grid steps
//   N     - number of time steps
//   rannacher_steps - number of fully-implicit startup steps
BSResult fd_price(double S, double K, double r, double q, double sigma,
                   double T, OptionType type, ExerciseStyle style, int M,
                   int N, int rannacher_steps = 2);

} // namespace optionspricer
