#include "optionspricer/binomial_tree.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace optionspricer {

double binomial_price(double S, double K, double r, double q, double sigma,
                       double T, OptionType type, ExerciseStyle style,
                       std::uint32_t steps) {
    const double dt = T / static_cast<double>(steps);
    const double u = std::exp(sigma * std::sqrt(dt));
    const double d = 1.0 / u;
    const double growth = std::exp((r - q) * dt);
    const double p = (growth - d) / (u - d);
    const double discount = std::exp(-r * dt);

    // Terminal asset prices and payoffs, indexed by number of up-moves.
    std::vector<double> values(steps + 1);
    for (std::uint32_t i = 0; i <= steps; ++i) {
        const double S_T = S * std::pow(u, static_cast<double>(i)) *
                            std::pow(d, static_cast<double>(steps - i));
        values[i] = (type == OptionType::Call) ? std::max(S_T - K, 0.0)
                                                 : std::max(K - S_T, 0.0);
    }

    // Backward induction.
    for (std::uint32_t step = steps; step-- > 0;) {
        for (std::uint32_t i = 0; i <= step; ++i) {
            const double continuation =
                discount * (p * values[i + 1] + (1.0 - p) * values[i]);

            if (style == ExerciseStyle::American) {
                const double S_node =
                    S * std::pow(u, static_cast<double>(i)) *
                    std::pow(d, static_cast<double>(step - i));
                const double exercise = (type == OptionType::Call)
                                             ? std::max(S_node - K, 0.0)
                                             : std::max(K - S_node, 0.0);
                values[i] = std::max(continuation, exercise);
            } else {
                values[i] = continuation;
            }
        }
    }

    return values[0];
}

} // namespace optionspricer
