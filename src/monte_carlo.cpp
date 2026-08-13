#include "optionspricer/monte_carlo.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace optionspricer {

namespace {

inline double payoff(double S_T, double K, OptionType type) {
    if (type == OptionType::Call) {
        return std::max(S_T - K, 0.0);
    }
    return std::max(K - S_T, 0.0);
}

// Sample mean and (unbiased) sample variance from running sums.
struct MeanVar {
    double mean = 0.0;
    double variance = 0.0;
};

MeanVar mean_and_variance(double sum, double sumsq, std::uint64_t n) {
    MeanVar mv;
    mv.mean = sum / static_cast<double>(n);
    if (n > 1) {
        const double sumsq_over_n = sumsq / static_cast<double>(n);
        mv.variance = (sumsq_over_n - mv.mean * mv.mean) *
                      static_cast<double>(n) / static_cast<double>(n - 1);
        mv.variance = std::max(mv.variance, 0.0);
    }
    return mv;
}

} // namespace

MCResult monte_carlo_price(double S, double K, double r, double q,
                            double sigma, double T, OptionType type,
                            std::uint64_t n_paths, bool antithetic,
                            std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double vol_sqrtT = sigma * std::sqrt(T);
    const double discount = std::exp(-r * T);

    // Pooled sums over every individual discounted payoff drawn, used to
    // estimate what plain (non-antithetic) MC variance would have been for
    // the same total number of normal draws.
    double sum_individual = 0.0;
    double sumsq_individual = 0.0;
    std::uint64_t n_individual = 0;

    // Sums over the actual estimator samples (paired averages when
    // antithetic, single draws otherwise).
    double sum_sample = 0.0;
    double sumsq_sample = 0.0;

    for (std::uint64_t i = 0; i < n_paths; ++i) {
        const double z = normal(rng);
        const double ST1 = S * std::exp(drift + vol_sqrtT * z);
        const double pay1 = discount * payoff(ST1, K, type);

        sum_individual += pay1;
        sumsq_individual += pay1 * pay1;
        ++n_individual;

        if (antithetic) {
            const double ST2 = S * std::exp(drift - vol_sqrtT * z);
            const double pay2 = discount * payoff(ST2, K, type);

            sum_individual += pay2;
            sumsq_individual += pay2 * pay2;
            ++n_individual;

            const double pair_avg = 0.5 * (pay1 + pay2);
            sum_sample += pair_avg;
            sumsq_sample += pair_avg * pair_avg;
        } else {
            sum_sample += pay1;
            sumsq_sample += pay1 * pay1;
        }
    }

    const MeanVar sample_mv = mean_and_variance(sum_sample, sumsq_sample, n_paths);

    MCResult res;
    res.price = sample_mv.mean;
    const double var_of_mean = sample_mv.variance / static_cast<double>(n_paths);
    res.std_error = std::sqrt(var_of_mean);
    res.ci_lower = res.price - 1.96 * res.std_error;
    res.ci_upper = res.price + 1.96 * res.std_error;

    if (antithetic) {
        const MeanVar indiv_mv =
            mean_and_variance(sum_individual, sumsq_individual, n_individual);
        const double naive_var_of_mean =
            indiv_mv.variance / static_cast<double>(n_individual);
        res.variance_reduction =
            (var_of_mean > 0.0) ? naive_var_of_mean / var_of_mean : 1.0;
    } else {
        res.variance_reduction = 1.0;
    }

    return res;
}

} // namespace optionspricer
