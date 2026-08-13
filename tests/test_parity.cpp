#include "optionspricer/black_scholes.hpp"
#include "test_utils.hpp"

#include <cmath>
#include <vector>

using namespace optionspricer;
using testutils::check_near;

int main() {
    const std::vector<double> moneyness = {0.8, 1.0, 1.2};
    const std::vector<double> vols = {0.15, 0.30};
    const std::vector<double> maturities = {0.25, 1.0, 2.0};
    const std::vector<double> rates = {0.03, -0.01};
    const double K = 100.0;
    const double q = 0.01;

    int n_checked = 0;
    for (double m : moneyness) {
        for (double sigma : vols) {
            for (double T : maturities) {
                for (double r : rates) {
                    const double S = m * K;
                    const double call =
                        black_scholes(S, K, r, q, sigma, T, OptionType::Call)
                            .price;
                    const double put =
                        black_scholes(S, K, r, q, sigma, T, OptionType::Put)
                            .price;
                    const double parity_rhs =
                        S * std::exp(-q * T) - K * std::exp(-r * T);
                    check_near(call - put, parity_rhs, 1e-6,
                               "put-call parity: S=" + std::to_string(S) +
                                   " sigma=" + std::to_string(sigma) +
                                   " T=" + std::to_string(T) +
                                   " r=" + std::to_string(r));
                    ++n_checked;
                }
            }
        }
    }

    check_near(static_cast<double>(n_checked), 36.0, 0.0,
               "grid covers the expected number of scenarios");

    return testutils::report("test_parity");
}
