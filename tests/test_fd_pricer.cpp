#include "optionspricer/black_scholes.hpp"
#include "optionspricer/fd_pricer.hpp"
#include "test_utils.hpp"

#include <cmath>
#include <cstdio>

using namespace optionspricer;
using testutils::check;
using testutils::check_near;

int main() {
    const double S0 = 100, K = 100, r = 0.05, q = 0.02, sigma = 0.20, T = 1.0;

    std::printf("European put, S0=%.0f K=%.0f r=%.2f q=%.2f sig=%.2f T=%.1f\n\n",
                S0, K, r, q, sigma, T);
    const double exact = black_scholes(S0, K, r, q, sigma, T, OptionType::Put).price;
    std::printf("  analytic Black-Scholes = %.10f\n\n", exact);

    // Convergence study: doubling M and N together should roughly quarter
    // the error, since Crank-Nicolson is second order in both space and
    // time (i.e. the abs-error ratio between successive refinements should
    // land close to 4.0).
    std::printf("  %6s %6s  %14s  %12s  %7s\n", "M", "N", "CN price", "abs error", "ratio");
    double prev = 0.0;
    for (int k = 0; k < 5; ++k) {
        const int M = 100 << k, N = 100 << k;
        const BSResult R = fd_price(S0, K, r, q, sigma, T, OptionType::Put,
                                     ExerciseStyle::European, M, N);
        const double err = std::fabs(R.price - exact);
        if (k == 0) {
            std::printf("  %6d %6d  %14.10f  %12.3e  %7s\n", M, N, R.price, err, "-");
        } else {
            const double ratio = prev / err;
            std::printf("  %6d %6d  %14.10f  %12.3e  %7.2f\n", M, N, R.price, err, ratio);
            check_near(ratio, 4.0, 0.5,
                       "CN abs-error ratio near 4.0 (M=N=" + std::to_string(M) + ")");
        }
        prev = err;
    }

    // Put-call parity on the grid: C - P = S*exp(-qT) - K*exp(-rT)
    const BSResult Pc = fd_price(S0, K, r, q, sigma, T, OptionType::Call,
                                  ExerciseStyle::European, 800, 800);
    const BSResult Pp = fd_price(S0, K, r, q, sigma, T, OptionType::Put,
                                  ExerciseStyle::European, 800, 800);
    const double parity_lhs = Pc.price - Pp.price;
    const double parity_rhs = S0 * std::exp(-q * T) - K * std::exp(-r * T);
    const double parity_residual = std::fabs(parity_lhs - parity_rhs);
    std::printf("\n  put-call parity residual (M=N=800) = %.3e\n", parity_residual);
    check(parity_residual < 1e-6,
          "put-call parity residual on the FD grid stays within 1e-6");

    // Greeks straight off the grid, no extra solves, checked against the
    // existing closed-form Black-Scholes implementation.
    const BSResult analytic_put =
        black_scholes(S0, K, r, q, sigma, T, OptionType::Put);
    std::printf("  delta: grid %.8f vs analytic %.8f  (err %.2e)\n",
                Pp.delta, analytic_put.delta, std::fabs(Pp.delta - analytic_put.delta));
    std::printf("  gamma: grid %.8f vs analytic %.8f  (err %.2e)\n",
                Pp.gamma, analytic_put.gamma, std::fabs(Pp.gamma - analytic_put.gamma));
    check_near(Pp.delta, analytic_put.delta, 1e-3, "FD delta matches analytic delta");
    check_near(Pp.gamma, analytic_put.gamma, 1e-3, "FD gamma matches analytic gamma");

    // American put: no closed form, so check the early-exercise premium is
    // positive relative to the FD European put on the same grid.
    const BSResult Am = fd_price(S0, K, r, q, sigma, T, OptionType::Put,
                                  ExerciseStyle::American, 800, 800);
    std::printf("\n  American put = %.8f   European put = %.8f   premium = %.8f\n",
                Am.price, Pp.price, Am.price - Pp.price);
    check(Am.price >= Pp.price - 1e-10,
          "American put is worth at least the European put on the FD grid");

    // Rannacher on/off, gamma at the strike, coarse grid: shows the ringing
    // that pure Crank-Nicolson produces against the kink in the payoff.
    const BSResult g_on = fd_price(S0, K, r, q, sigma, T, OptionType::Put,
                                    ExerciseStyle::European, 200, 25, 2);
    const BSResult g_off = fd_price(S0, K, r, q, sigma, T, OptionType::Put,
                                     ExerciseStyle::European, 200, 25, 0);
    std::printf("\n  gamma at K (M=200,N=25): Rannacher %.6f | pure CN %.6f | analytic %.6f\n",
                g_on.gamma, g_off.gamma, analytic_put.gamma);

    return testutils::report("test_fd_pricer");
}
