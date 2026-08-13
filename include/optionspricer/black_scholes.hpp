#pragma once

#include "optionspricer/types.hpp"

namespace optionspricer {

// Price plus the five standard Greeks, all in closed form.
//
// Conventions:
//   - vega is per unit change in volatility (i.e. per 1.00, not per 1%).
//   - rho is per unit change in the risk-free rate (i.e. per 1.00, not per 1%).
//   - theta is per year (not per calendar day).
struct BSResult {
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double vega = 0.0;
    double theta = 0.0;
    double rho = 0.0;
};

// European option price and Greeks under Black-Scholes-Merton.
//   S     - spot price
//   K     - strike price
//   r     - continuously compounded risk-free rate (may be negative)
//   q     - continuously compounded dividend yield
//   sigma - volatility (annualized)
//   T     - time to maturity in years
BSResult black_scholes(double S, double K, double r, double q, double sigma,
                        double T, OptionType type);

} // namespace optionspricer
