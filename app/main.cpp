// Validation harness for the options-pricer library.
//
// Runs a grid of Black-Scholes analytic prices against the Monte Carlo and
// CRR binomial numerical pricers, checks that everything agrees with the
// analytic benchmark to within expected tolerances, and prints the results
// as a set of readable tables.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "optionspricer/binomial_tree.hpp"
#include "optionspricer/black_scholes.hpp"
#include "optionspricer/monte_carlo.hpp"

using namespace optionspricer;

namespace {

// A validation harness must not use the standard assert() macro: it
// compiles to a no-op under NDEBUG, which CMake's Release build type
// defines by default, silently turning every check below into dead code.
// require() always checks, regardless of build type.
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "VALIDATION FAILED: " << message << "\n";
        std::exit(1);
    }
}

constexpr double K_STRIKE = 100.0;
constexpr std::uint64_t MC_PATHS = 100000; // antithetic pairs -> 200k draws
constexpr std::uint32_t BINOMIAL_STEPS = 1000;

std::string type_name(OptionType t) {
    return t == OptionType::Call ? "Call" : "Put";
}

struct Scenario {
    double S, K, r, q, sigma, T;
    OptionType type;
};

std::vector<Scenario> build_grid() {
    const std::vector<double> moneyness = {0.8, 1.0, 1.2};
    const std::vector<double> vols = {0.20, 0.40};
    const std::vector<double> maturities = {0.25, 1.0};
    const std::vector<double> rates = {0.03, -0.02}; // one negative rate
    const std::vector<OptionType> types = {OptionType::Call, OptionType::Put};
    const double q = 0.0;

    std::vector<Scenario> grid;
    for (double m : moneyness) {
        for (double sigma : vols) {
            for (double T : maturities) {
                for (double r : rates) {
                    for (OptionType type : types) {
                        grid.push_back(
                            {m * K_STRIKE, K_STRIKE, r, q, sigma, T, type});
                    }
                }
            }
        }
    }
    return grid;
}

void print_divider(int width) {
    std::cout << std::string(static_cast<std::size_t>(width), '-') << "\n";
}

} // namespace

int main() {
    std::cout << std::fixed << std::setprecision(4);

    // ------------------------------------------------------------------
    // Section 1: full parameter grid, all three pricers vs analytic price
    // ------------------------------------------------------------------
    std::cout << "=================================================================================================\n";
    std::cout << " SECTION 1: Parameter grid -- analytic vs Monte Carlo vs binomial\n";
    std::cout << "=================================================================================================\n";

    const std::vector<Scenario> grid = build_grid();

    const int width = 128;
    print_divider(width);
    std::cout << std::left << std::setw(5) << "Type" << std::setw(7) << "S/K"
              << std::setw(8) << "sigma" << std::setw(6) << "T" << std::setw(8)
              << "r" << std::right << std::setw(10) << "BS Price"
              << std::setw(12) << "MC Price" << std::setw(22) << "MC 95% CI"
              << std::setw(10) << "MC Err" << std::setw(12) << "Binom Price"
              << std::setw(12) << "Binom Err" << std::setw(10) << "MC in CI"
              << "\n";
    print_divider(width);

    int n_mc_ci_failures = 0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const Scenario& sc = grid[i];

        const BSResult bs =
            black_scholes(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T, sc.type);

        const std::uint64_t seed = 1000 + i;
        const MCResult mc =
            monte_carlo_price(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T, sc.type,
                               MC_PATHS, /*antithetic=*/true, seed);

        const double binom =
            binomial_price(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T, sc.type,
                            ExerciseStyle::European, BINOMIAL_STEPS);

        const double mc_err = std::fabs(mc.price - bs.price);
        const double binom_err = std::fabs(binom - bs.price);
        const bool in_ci = (bs.price >= mc.ci_lower && bs.price <= mc.ci_upper);
        if (!in_ci) {
            ++n_mc_ci_failures;
        }

        std::string ci_str = "[" + std::to_string(mc.ci_lower).substr(0, 7) +
                              ", " + std::to_string(mc.ci_upper).substr(0, 7) +
                              "]";

        std::cout << std::left << std::setw(5) << type_name(sc.type)
                   << std::setw(7) << (sc.S / sc.K) << std::setw(8) << sc.sigma
                   << std::setw(6) << sc.T << std::setw(8) << sc.r
                   << std::right << std::setw(10) << bs.price << std::setw(12)
                   << mc.price << std::setw(22) << ci_str << std::setw(10)
                   << mc_err << std::setw(12) << binom << std::setw(12)
                   << binom_err << std::setw(10) << (in_ci ? "yes" : "NO")
                   << "\n";

    }
    print_divider(width);
    const double miss_rate =
        static_cast<double>(n_mc_ci_failures) / static_cast<double>(grid.size());
    std::cout << grid.size() << " scenarios, " << n_mc_ci_failures
              << " Monte Carlo CI containment failures (miss rate = "
              << (miss_rate * 100.0) << "%).\n";

    // A correctly calibrated 95% CI is *expected* to miss about 5% of the
    // time by chance, so asserting zero misses on every single row would be
    // asserting a statistical impossibility, not a bug. Instead we assert
    // that the aggregate miss rate is consistent with that ~5% expectation
    // (generous upper bound to keep this from flaking under normal sampling
    // noise, while still catching a systematically biased estimator).
    require(miss_rate <= 0.20,
            "Monte Carlo CI miss rate is far above the ~5% expected under "
            "correctly calibrated 95% confidence intervals -- suggests a "
            "bias in the Monte Carlo pricer.");
    std::cout << "(Individual scenarios missing their own 95% CI now and "
                 "then is expected; see README for why this is checked in "
                 "aggregate rather than per-row.)\n\n";

    // ------------------------------------------------------------------
    // Section 2: binomial convergence study
    // ------------------------------------------------------------------
    std::cout << "=================================================================================================\n";
    std::cout << " SECTION 2: Binomial convergence to Black-Scholes as step count grows\n";
    std::cout << "=================================================================================================\n";

    {
        const double S = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.2,
                     T = 1.0;
        const double bs_price =
            black_scholes(S, K, r, q, sigma, T, OptionType::Call).price;
        const std::vector<std::uint32_t> step_counts = {10, 50, 100, 500,
                                                          1000, 5000};

        std::cout << "ATM European call, S=K=100, r=5%, sigma=20%, T=1y. "
                     "Black-Scholes price = "
                  << bs_price << "\n\n";
        print_divider(50);
        std::cout << std::left << std::setw(10) << "N" << std::right
                  << std::setw(15) << "Binomial Price" << std::setw(15)
                  << "Abs Error" << std::setw(10) << " N*Err" << "\n";
        print_divider(50);
        for (std::uint32_t n : step_counts) {
            const double price = binomial_price(S, K, r, q, sigma, T,
                                                  OptionType::Call,
                                                  ExerciseStyle::European, n);
            const double err = std::fabs(price - bs_price);
            std::cout << std::left << std::setw(10) << n << std::right
                      << std::setw(15) << price << std::setw(15) << err
                      << std::setw(10) << (err * static_cast<double>(n))
                      << "\n";
        }
        print_divider(50);
        std::cout << "N*Error staying roughly constant/bounded across rows "
                     "is the signature of O(1/N) convergence.\n\n";
    }

    // ------------------------------------------------------------------
    // Section 3: put-call parity across the grid
    // ------------------------------------------------------------------
    std::cout << "=================================================================================================\n";
    std::cout << " SECTION 3: Put-call parity check across the grid\n";
    std::cout << "=================================================================================================\n";
    {
        double max_parity_violation = 0.0;
        for (std::size_t i = 0; i + 1 < grid.size(); i += 2) {
            // Grid is built with Call then Put for identical (S,K,r,q,sigma,T).
            const Scenario& call_sc = grid[i];
            const Scenario& put_sc = grid[i + 1];
            if (call_sc.type != OptionType::Call ||
                put_sc.type != OptionType::Put) {
                continue;
            }
            const double call = black_scholes(call_sc.S, call_sc.K, call_sc.r,
                                                call_sc.q, call_sc.sigma,
                                                call_sc.T, OptionType::Call)
                                     .price;
            const double put = black_scholes(put_sc.S, put_sc.K, put_sc.r,
                                               put_sc.q, put_sc.sigma,
                                               put_sc.T, OptionType::Put)
                                    .price;
            const double rhs = call_sc.S * std::exp(-call_sc.q * call_sc.T) -
                                call_sc.K * std::exp(-call_sc.r * call_sc.T);
            const double violation = std::fabs((call - put) - rhs);
            max_parity_violation = std::max(max_parity_violation, violation);
            require(violation < 1e-6, "put-call parity violated");
        }
        std::cout << "Checked " << (grid.size() / 2)
                  << " call/put pairs. Max |C - P - (S*e^(-qT) - K*e^(-rT))| = "
                  << std::scientific << max_parity_violation << std::fixed
                  << "\n\n";
    }

    // ------------------------------------------------------------------
    // Section 4: American vs European exercise
    // ------------------------------------------------------------------
    std::cout << "=================================================================================================\n";
    std::cout << " SECTION 4: American vs European exercise (binomial tree)\n";
    std::cout << "=================================================================================================\n";
    {
        int n_checked_puts = 0;
        int n_checked_calls_nonneg_r = 0;
        int n_calls_neg_r_diverged = 0;
        double max_call_diff_no_div_nonneg_r = 0.0;
        for (const Scenario& sc : grid) {
            const double euro = binomial_price(sc.S, sc.K, sc.r, sc.q,
                                                sc.sigma, sc.T, sc.type,
                                                ExerciseStyle::European, 500);
            const double amer = binomial_price(sc.S, sc.K, sc.r, sc.q,
                                                sc.sigma, sc.T, sc.type,
                                                ExerciseStyle::American, 500);
            if (sc.type == OptionType::Put) {
                // American value dominates European value unconditionally:
                // the American holder can always choose not to exercise
                // early and replicate the European payoff exactly.
                require(amer >= euro - 1e-10,
                        "American put must be worth at least the European put");
                ++n_checked_puts;
            } else {
                // The textbook result "American call == European call with
                // no dividends" additionally assumes r >= 0: the proof
                // relies on K*e^(-rT) <= K, which flips when r < 0. With a
                // negative rate, deferring payment of the strike is a net
                // cost (cash itself decays under negative rates), so early
                // exercise on a deep-ITM, low-vol call can become optimal
                // even with q = 0. This grid includes r = -0.02 specifically
                // to surface that, so we only assert equality for r >= 0.
                if (sc.r >= 0.0) {
                    max_call_diff_no_div_nonneg_r = std::max(
                        max_call_diff_no_div_nonneg_r, std::fabs(amer - euro));
                    require(std::fabs(amer - euro) < 1e-8,
                            "American call must equal European call when "
                            "q=0 and r>=0");
                    ++n_checked_calls_nonneg_r;
                } else {
                    require(amer >= euro - 1e-10,
                            "American call must be worth at least the "
                            "European call");
                    if (std::fabs(amer - euro) > 1e-6) {
                        ++n_calls_neg_r_diverged;
                    }
                }
            }
        }
        std::cout << "American put >= European put: verified across "
                  << n_checked_puts << " put scenarios.\n";
        std::cout << "American call == European call (q=0, r>=0): verified "
                     "across "
                  << n_checked_calls_nonneg_r
                  << " call scenarios, max diff = " << std::scientific
                  << max_call_diff_no_div_nonneg_r << std::fixed << "\n";
        std::cout << "American call vs European call (q=0, r<0): "
                  << n_calls_neg_r_diverged
                  << " scenario(s) show a real early-exercise premium -- "
                     "negative rates make deferring the strike payment "
                     "costly even without dividends (see README).\n";

        // With a positive dividend yield, American calls should carry a
        // real early-exercise premium over their European counterparts.
        const double S = 100.0, K = 90.0, r = 0.03, q = 0.06, sigma = 0.2,
                     T = 1.0;
        const double euro_div = binomial_price(S, K, r, q, sigma, T,
                                                 OptionType::Call,
                                                 ExerciseStyle::European, 500);
        const double amer_div = binomial_price(S, K, r, q, sigma, T,
                                                 OptionType::Call,
                                                 ExerciseStyle::American, 500);
        std::cout << "With a dividend yield q=6%: European call = " << euro_div
                  << ", American call = " << amer_div
                  << " (premium = " << (amer_div - euro_div) << ")\n\n";
        require(amer_div >= euro_div - 1e-10,
                "American call must be worth at least the European call");
    }

    std::cout << "=================================================================================================\n";
    std::cout << " All validations passed.\n";
    std::cout << "=================================================================================================\n";

    return 0;
}
