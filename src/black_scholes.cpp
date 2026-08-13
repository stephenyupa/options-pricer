#include "optionspricer/black_scholes.hpp"
#include "optionspricer/normal_dist.hpp"

#include <cmath>

namespace optionspricer {

BSResult black_scholes(double S, double K, double r, double q, double sigma,
                        double T, OptionType type) {
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) /
                       (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;

    const double disc_r = std::exp(-r * T);
    const double disc_q = std::exp(-q * T);
    const double pdf_d1 = norm_pdf(d1);

    BSResult res;

    if (type == OptionType::Call) {
        const double Nd1 = norm_cdf(d1);
        const double Nd2 = norm_cdf(d2);
        res.price = S * disc_q * Nd1 - K * disc_r * Nd2;
        res.delta = disc_q * Nd1;
        res.theta = -(S * disc_q * pdf_d1 * sigma) / (2.0 * sqrtT) -
                    r * K * disc_r * Nd2 + q * S * disc_q * Nd1;
        res.rho = K * T * disc_r * Nd2;
    } else {
        const double Nnd1 = norm_cdf(-d1);
        const double Nnd2 = norm_cdf(-d2);
        res.price = K * disc_r * Nnd2 - S * disc_q * Nnd1;
        res.delta = disc_q * (norm_cdf(d1) - 1.0);
        res.theta = -(S * disc_q * pdf_d1 * sigma) / (2.0 * sqrtT) +
                    r * K * disc_r * Nnd2 - q * S * disc_q * Nnd1;
        res.rho = -K * T * disc_r * Nnd2;
    }

    // gamma and vega are identical for calls and puts.
    res.gamma = disc_q * pdf_d1 / (S * sigma * sqrtT);
    res.vega = S * disc_q * pdf_d1 * sqrtT;

    return res;
}

} // namespace optionspricer
