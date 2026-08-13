#pragma once

#include <cstdint>

#include "optionspricer/types.hpp"

namespace optionspricer {

// Cox-Ross-Rubinstein binomial tree price for a European or American option.
//   S     - spot price
//   K     - strike price
//   r     - continuously compounded risk-free rate (may be negative)
//   q     - continuously compounded dividend yield
//   sigma - volatility (annualized)
//   T     - time to maturity in years
//   steps - number of time steps in the tree
double binomial_price(double S, double K, double r, double q, double sigma,
                       double T, OptionType type, ExerciseStyle style,
                       std::uint32_t steps);

} // namespace optionspricer
