// C-linkage shim over optionspricer::binomial_price and black_scholes,
// built into a shared library so the Julia test suite can ccall the *exact*
// C++ implementations directly, rather than re-deriving expected values by
// hand. The Black-Scholes export exists only as a convergence benchmark for
// the binomial test grid.
#include "optionspricer/binomial_tree.hpp"
#include "optionspricer/black_scholes.hpp"

using namespace optionspricer;

extern "C" double binomial_price_c(double S, double K, double r, double q,
                                    double sigma, double T, int type,
                                    int style, unsigned int steps) {
    return binomial_price(S, K, r, q, sigma, T,
                           type == 0 ? OptionType::Call : OptionType::Put,
                           style == 0 ? ExerciseStyle::European
                                      : ExerciseStyle::American,
                           steps);
}

extern "C" double black_scholes_price_c(double S, double K, double r,
                                         double q, double sigma, double T,
                                         int type) {
    return black_scholes(S, K, r, q, sigma, T,
                          type == 0 ? OptionType::Call : OptionType::Put)
        .price;
}
