#pragma once

namespace optionspricer {

// Standard normal probability density function.
double norm_pdf(double x);

// Standard normal cumulative distribution function, computed via
// N(x) = 0.5 * erfc(-x / sqrt(2)) for full double precision accuracy.
double norm_cdf(double x);

} // namespace optionspricer
