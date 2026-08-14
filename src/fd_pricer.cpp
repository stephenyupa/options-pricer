#include "optionspricer/fd_pricer.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace optionspricer {

namespace {

// Thomas algorithm: solve a tridiagonal system in O(n). `lo`, `di`, `up` are
// the sub-, main-, and super-diagonals; `d` is the RHS. `di` is copied
// because the sweep destroys it.
std::vector<double> thomas(const std::vector<double>& lo, std::vector<double> di,
                            const std::vector<double>& up,
                            std::vector<double> d) {
    const int n = (int)di.size();
    for (int i = 1; i < n; ++i) {
        const double w = lo[i] / di[i - 1];
        di[i] -= w * up[i - 1];
        d[i]  -= w * d[i - 1];
    }
    std::vector<double> x(n);
    x[n - 1] = d[n - 1] / di[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = (d[i] - up[i] * x[i + 1]) / di[i];
    return x;
}

// Projected SOR for the linear complementarity problem that American exercise
// induces: solve the tridiagonal system subject to V >= payoff at every node.
std::vector<double> psor(const std::vector<double>& lo, const std::vector<double>& di,
                         const std::vector<double>& up, const std::vector<double>& d,
                         const std::vector<double>& payoff, std::vector<double> x,
                         double omega = 1.2, double tol = 1e-12, int max_it = 20000) {
    const int n = (int)di.size();
    for (int it = 0; it < max_it; ++it) {
        double err = 0.0;
        for (int i = 0; i < n; ++i) {
            double resid = d[i] - di[i] * x[i];
            if (i > 0)     resid -= lo[i] * x[i - 1];
            if (i < n - 1) resid -= up[i] * x[i + 1];
            const double cand = std::max(payoff[i], x[i] + omega * resid / di[i]);
            err = std::max(err, std::fabs(cand - x[i]));
            x[i] = cand;
        }
        if (err < tol) break;
    }
    return x;
}

} // namespace

BSResult fd_price(double S0, double K, double r, double q, double sigma,
                   double T, OptionType type, ExerciseStyle style, int M,
                   int N, int rannacher_steps) {
    const bool is_call = (type == OptionType::Call);
    const bool american = (style == ExerciseStyle::American);

    // Domain: wide enough that the far boundary is not felt at S0. Five
    // standard deviations of terminal log-return plus drift is comfortable.
    double S_max = std::max(3.0 * K, S0 * std::exp((r - q) * T + 5.0 * sigma * std::sqrt(T)));

    // Pin S0 to a node so the reported price needs no interpolation, and so
    // delta/gamma come from clean central differences.
    int i0 = std::max(1, (int)std::llround(S0 / (S_max / M)));
    const double dS = S0 / i0;
    S_max = M * dS;
    if (i0 >= M - 1) { i0 = M / 2; }  // degenerate guard

    const double dt = T / N;

    std::vector<double> S(M + 1), V(M + 1);
    for (int i = 0; i <= M; ++i) {
        S[i] = i * dS;
        V[i] = std::max(is_call ? S[i] - K : K - S[i], 0.0);
    }

    // Spatial operator coefficients (independent of dS, as derived above).
    std::vector<double> a(M), b(M), c(M);
    for (int i = 1; i < M; ++i) {
        const double i2 = (double)i * i;
        a[i] = 0.5 * sigma * sigma * i2 - 0.5 * (r - q) * i;
        b[i] = -sigma * sigma * i2 - r;
        c[i] = 0.5 * sigma * sigma * i2 + 0.5 * (r - q) * i;
    }

    // Intrinsic value, used as the American exercise constraint.
    std::vector<double> intrinsic(M - 1);
    for (int i = 1; i < M; ++i)
        intrinsic[i - 1] = std::max(is_call ? S[i] - K : K - S[i], 0.0);

    std::vector<double> lo(M - 1), di(M - 1), up(M - 1), rhs(M - 1);

    for (int n = 0; n < N; ++n) {
        const double tau_new = (n + 1) * dt;
        // Rannacher: fully implicit for the first few steps, CN thereafter.
        const double th = (n < rannacher_steps) ? 1.0 : 0.5;

        // Dirichlet boundaries at the new time level.
        double V0, VM;
        if (is_call) {
            V0 = 0.0;
            VM = american ? std::max(S_max - K, S_max * std::exp(-q * tau_new) - K * std::exp(-r * tau_new))
                          : S_max * std::exp(-q * tau_new) - K * std::exp(-r * tau_new);
        } else {
            V0 = american ? K : K * std::exp(-r * tau_new);
            VM = 0.0;
        }

        for (int i = 1; i < M; ++i) {
            const int k = i - 1;
            lo[k] = -th * dt * a[i];
            di[k] = 1.0 - th * dt * b[i];
            up[k] = -th * dt * c[i];

            const double ex = (1.0 - th) * dt;
            rhs[k] = ex * a[i] * V[i - 1] + (1.0 + ex * b[i]) * V[i] + ex * c[i] * V[i + 1];
        }
        // Move the known boundary terms of the implicit operator to the RHS.
        rhs[0]     += th * dt * a[1]     * V0;
        rhs[M - 2] += th * dt * c[M - 1] * VM;

        std::vector<double> guess(V.begin() + 1, V.begin() + M);
        std::vector<double> sol = american
            ? psor(lo, di, up, rhs, intrinsic, guess)
            : thomas(lo, di, up, rhs);

        V[0] = V0;
        V[M] = VM;
        for (int i = 1; i < M; ++i) V[i] = sol[i - 1];
    }

    BSResult out;
    out.price = V[i0];
    out.delta = (V[i0 + 1] - V[i0 - 1]) / (2.0 * dS);
    out.gamma = (V[i0 + 1] - 2.0 * V[i0] + V[i0 - 1]) / (dS * dS);
    return out;
}

} // namespace optionspricer
