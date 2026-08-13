#include "optionspricer/normal_dist.hpp"

#include <cmath>

namespace optionspricer {

double norm_pdf(double x) {
    static const double inv_sqrt_2pi = 0.3989422804014327; // 1 / sqrt(2*pi)
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

} // namespace optionspricer
