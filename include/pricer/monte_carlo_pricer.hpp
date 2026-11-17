#pragma once

#include <cstddef>
#include <cstdint>

#include "pricer/pricer.hpp"

namespace pricer {

/// Variance reduction technique used by MonteCarloPricer.
enum class VarianceReduction {
    /// Plain Monte Carlo: average of N independent discounted payoffs.
    None,
    /// Antithetic variates: N/2 pairs (Z, -Z), each pair averaged. Same number of payoff
    /// evaluations as plain MC, so results are compared at equal cost.
    Antithetic,
    /// Control variate X = e^{-rT} S_T with known mean S_0 e^{-qT}:
    /// estimator mean(h) - beta (mean(X) - S_0 e^{-qT}), beta = Cov(h, X) / Var(X)
    /// estimated on the same paths.
    ControlVariate
};

/// Monte Carlo pricer for European options under Black-Scholes (geometric Brownian motion).
///
/// Uses the exact solution of the SDE, so a single time step is enough and there is
/// no discretisation bias:
///
///     S_T = S_0 exp((r - q - sigma^2/2) T + sigma sqrt(T) Z),   Z ~ N(0, 1)
///
/// The estimator is the discounted sample mean of the payoff (optionally with a variance
/// reduction technique). The standard error is s / sqrt(n) with s the sample standard
/// deviation of the n i.i.d. samples actually averaged (pairs for antithetic variates),
/// and the 95% confidence interval is price +/- 1.96 SE (central limit theorem).
///
/// Works for any Payoff. Pricing is const and deterministic: each call re-seeds a fresh
/// generator with the stored seed, so repeated calls (for instance with bumped market
/// data) reuse the same random numbers.
class MonteCarloPricer final : public Pricer {
public:
    static constexpr std::uint64_t kDefaultSeed = 42;

    /// @param numPaths  number of simulated terminal prices (payoff evaluations); at least 2,
    ///                  and at least 4 for antithetic variates (2 pairs)
    /// @param seed      seed of the random number generator
    /// @param reduction variance reduction technique
    /// @throws std::invalid_argument if numPaths is too small.
    explicit MonteCarloPricer(std::size_t numPaths, std::uint64_t seed = kDefaultSeed,
                              VarianceReduction reduction = VarianceReduction::None);

    /// @throws std::invalid_argument for American exercise.
    [[nodiscard]] PricingResult price(const Option& option,
                                      const MarketData& market) const override;

    [[nodiscard]] std::size_t numPaths() const noexcept { return numPaths_; }
    [[nodiscard]] std::uint64_t seed() const noexcept { return seed_; }
    [[nodiscard]] VarianceReduction reduction() const noexcept { return reduction_; }

    /// Two-sided 97.5% quantile of N(0,1), used for the 95% confidence interval.
    static constexpr double kZ95 = 1.959963984540054;

private:
    std::size_t numPaths_;
    std::uint64_t seed_;
    VarianceReduction reduction_;
};

}  // namespace pricer
