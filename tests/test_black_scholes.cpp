#include "optionspricer/black_scholes.hpp"
#include "test_utils.hpp"

#include <cmath>

using namespace optionspricer;
using testutils::check_near;

namespace {

double price_only(double S, double K, double r, double q, double sigma,
                   double T, OptionType type) {
    return black_scholes(S, K, r, q, sigma, T, type).price;
}

// Checks the closed-form Greeks against central finite differences of the
// closed-form price itself, for both calls and puts.
void check_greeks_vs_finite_difference(double S, double K, double r, double q,
                                        double sigma, double T,
                                        OptionType type, double tol) {
    const BSResult res = black_scholes(S, K, r, q, sigma, T, type);
    const std::string tag = (type == OptionType::Call) ? "call" : "put";

    const double hS = 1e-3 * S;
    const double delta_fd = (price_only(S + hS, K, r, q, sigma, T, type) -
                              price_only(S - hS, K, r, q, sigma, T, type)) /
                             (2.0 * hS);
    check_near(res.delta, delta_fd, tol, tag + " delta vs finite difference");

    const double gamma_fd = (price_only(S + hS, K, r, q, sigma, T, type) -
                              2.0 * price_only(S, K, r, q, sigma, T, type) +
                              price_only(S - hS, K, r, q, sigma, T, type)) /
                             (hS * hS);
    check_near(res.gamma, gamma_fd, tol, tag + " gamma vs finite difference");

    const double hSig = 1e-4;
    const double vega_fd =
        (price_only(S, K, r, q, sigma + hSig, T, type) -
         price_only(S, K, r, q, sigma - hSig, T, type)) /
        (2.0 * hSig);
    check_near(res.vega, vega_fd, tol, tag + " vega vs finite difference");

    const double hT = 1e-5;
    const double theta_fd = -(price_only(S, K, r, q, sigma, T + hT, type) -
                               price_only(S, K, r, q, sigma, T - hT, type)) /
                             (2.0 * hT);
    check_near(res.theta, theta_fd, tol, tag + " theta vs finite difference");

    const double hR = 1e-5;
    const double rho_fd = (price_only(S, K, r + hR, q, sigma, T, type) -
                            price_only(S, K, r - hR, q, sigma, T, type)) /
                           (2.0 * hR);
    check_near(res.rho, rho_fd, tol, tag + " rho vs finite difference");
}

} // namespace

int main() {
    // Textbook reference: Hull's example, S=K=100, r=5%, sigma=20%, T=1y.
    check_near(price_only(100.0, 100.0, 0.05, 0.0, 0.20, 1.0, OptionType::Call),
               10.4506, 0.01, "ATM call price matches known reference");
    check_near(price_only(100.0, 100.0, 0.05, 0.0, 0.20, 1.0, OptionType::Put),
               5.5735, 0.01, "ATM put price matches known reference");

    // Put-call parity: C - P = S*exp(-qT) - K*exp(-rT).
    {
        const double S = 105.0, K = 95.0, r = 0.03, q = 0.01, sigma = 0.25,
                     T = 0.75;
        const double call = price_only(S, K, r, q, sigma, T, OptionType::Call);
        const double put = price_only(S, K, r, q, sigma, T, OptionType::Put);
        const double parity_rhs = S * std::exp(-q * T) - K * std::exp(-r * T);
        check_near(call - put, parity_rhs, 1e-8, "put-call parity holds");
    }

    // Greeks vs finite differences, across a few scenarios including a
    // negative rate and a nonzero dividend yield.
    check_greeks_vs_finite_difference(100.0, 100.0, 0.05, 0.0, 0.20, 1.0,
                                       OptionType::Call, 1e-3);
    check_greeks_vs_finite_difference(100.0, 100.0, 0.05, 0.0, 0.20, 1.0,
                                       OptionType::Put, 1e-3);
    check_greeks_vs_finite_difference(80.0, 100.0, -0.01, 0.02, 0.35, 0.25,
                                       OptionType::Call, 1e-3);
    check_greeks_vs_finite_difference(120.0, 100.0, -0.01, 0.02, 0.35, 2.0,
                                       OptionType::Put, 1e-3);

    return testutils::report("test_black_scholes");
}
