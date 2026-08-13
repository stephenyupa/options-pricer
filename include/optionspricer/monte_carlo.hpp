#pragma once

#include <cstdint>

#include "optionspricer/types.hpp"

namespace optionspricer {

struct MCResult {
    double price = 0.0;
    double std_error = 0.0;
    double ci_lower = 0.0; // 95% confidence interval, lower bound
    double ci_upper = 0.0; // 95% confidence interval, upper bound

    // Ratio of the naive-estimator variance to the achieved estimator
    // variance, both measured using the same total number of normal draws.
    // Only meaningful when antithetic variates were requested; left at 1.0
    // otherwise. > 1 means variance was reduced.
    double variance_reduction = 1.0;
};

// Monte Carlo price of a European option under geometric Brownian motion,
// simulated via the exact GBM solution:
//   S_T = S_0 * exp((r - q - sigma^2/2) * T + sigma * sqrt(T) * Z)
//
// n_paths is the number of independent samples averaged into the estimator.
// When antithetic is true, n_paths pairs (Z, -Z) are drawn (2*n_paths total
// normal draws) and each pair's averaged payoff counts as one sample, which
// is the standard antithetic-variates construction.
MCResult monte_carlo_price(double S, double K, double r, double q,
                            double sigma, double T, OptionType type,
                            std::uint64_t n_paths, bool antithetic,
                            std::uint64_t seed);

} // namespace optionspricer
