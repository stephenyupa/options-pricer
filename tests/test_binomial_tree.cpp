#include "optionspricer/binomial_tree.hpp"
#include "optionspricer/black_scholes.hpp"
#include "test_utils.hpp"

using namespace optionspricer;
using testutils::check;
using testutils::check_near;

int main() {
    // European binomial price converges to the Black-Scholes analytic price
    // as step count grows.
    {
        const double S = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.2,
                     T = 1.0;
        const double bs = black_scholes(S, K, r, q, sigma, T, OptionType::Call).price;
        const double tree_coarse =
            binomial_price(S, K, r, q, sigma, T, OptionType::Call,
                            ExerciseStyle::European, 25);
        const double tree_fine =
            binomial_price(S, K, r, q, sigma, T, OptionType::Call,
                            ExerciseStyle::European, 2000);

        check(std::fabs(tree_fine - bs) < std::fabs(tree_coarse - bs),
              "finer tree is closer to Black-Scholes than coarser tree");
        check_near(tree_fine, bs, 0.01,
                   "European binomial converges to Black-Scholes analytic price");
    }

    // American call == European call when there is no dividend (early
    // exercise is never optimal on a non-dividend-paying underlying).
    {
        const double S = 100.0, K = 90.0, r = 0.05, q = 0.0, sigma = 0.3,
                     T = 1.0;
        const double euro = binomial_price(S, K, r, q, sigma, T, OptionType::Call,
                                            ExerciseStyle::European, 500);
        const double amer = binomial_price(S, K, r, q, sigma, T, OptionType::Call,
                                            ExerciseStyle::American, 500);
        check_near(amer, euro, 1e-8,
                   "American call equals European call with zero dividends");
    }

    // American put is worth at least as much as the European put.
    {
        const double S = 100.0, K = 110.0, r = 0.04, q = 0.0, sigma = 0.25,
                     T = 1.0;
        const double euro = binomial_price(S, K, r, q, sigma, T, OptionType::Put,
                                            ExerciseStyle::European, 500);
        const double amer = binomial_price(S, K, r, q, sigma, T, OptionType::Put,
                                            ExerciseStyle::American, 500);
        check(amer >= euro - 1e-10,
              "American put is worth at least the European put");
        check(amer > euro + 1e-6,
              "American put early-exercise premium is strictly positive ITM");
    }

    // With a dividend yield, the American call can also carry an early
    // exercise premium over the European call.
    {
        const double S = 100.0, K = 90.0, r = 0.03, q = 0.06, sigma = 0.2,
                     T = 1.0;
        const double euro = binomial_price(S, K, r, q, sigma, T, OptionType::Call,
                                            ExerciseStyle::European, 500);
        const double amer = binomial_price(S, K, r, q, sigma, T, OptionType::Call,
                                            ExerciseStyle::American, 500);
        check(amer >= euro - 1e-10,
              "American call is worth at least the European call with dividends");
    }

    return testutils::report("test_binomial_tree");
}
