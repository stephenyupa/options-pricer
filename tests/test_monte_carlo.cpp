#include "optionspricer/black_scholes.hpp"
#include "optionspricer/monte_carlo.hpp"
#include "test_utils.hpp"

using namespace optionspricer;
using testutils::check;

int main() {
    const double S = 100.0, K = 105.0, r = 0.04, q = 0.0, sigma = 0.25,
                 T = 0.5;
    const double bs_price =
        black_scholes(S, K, r, q, sigma, T, OptionType::Call).price;

    // Plain Monte Carlo: the analytic price should fall inside the reported
    // 95% confidence interval.
    {
        const MCResult mc = monte_carlo_price(S, K, r, q, sigma, T,
                                               OptionType::Call, 300000,
                                               /*antithetic=*/false, 42);
        check(mc.ci_lower <= bs_price && bs_price <= mc.ci_upper,
              "plain MC confidence interval contains the analytic price");
        check(mc.ci_upper > mc.ci_lower, "confidence interval is non-degenerate");
    }

    // Antithetic Monte Carlo: CI should still contain the analytic price,
    // and variance reduction should be reported as an improvement (>1).
    {
        const MCResult mc = monte_carlo_price(S, K, r, q, sigma, T,
                                               OptionType::Call, 150000,
                                               /*antithetic=*/true, 42);
        check(mc.ci_lower <= bs_price && bs_price <= mc.ci_upper,
              "antithetic MC confidence interval contains the analytic price");
        check(mc.variance_reduction > 1.0,
              "antithetic variates reduce estimator variance");
    }

    // More paths should shrink the standard error.
    {
        const MCResult mc_small = monte_carlo_price(
            S, K, r, q, sigma, T, OptionType::Call, 5000, false, 7);
        const MCResult mc_large = monte_carlo_price(
            S, K, r, q, sigma, T, OptionType::Call, 500000, false, 7);
        check(mc_large.std_error < mc_small.std_error,
              "standard error shrinks as path count grows");
    }

    return testutils::report("test_monte_carlo");
}
