#pragma once

#include <cstddef>
#include <cstdint>

#include "pricer/pricer.hpp"

namespace pricer {

/// Monte Carlo pricer for European options under Black-Scholes (geometric Brownian motion).
///
/// Uses the exact solution of the SDE, so a single time step is enough and there is
/// no discretisation bias:
///
///     S_T = S_0 exp((r - q - sigma^2/2) T + sigma sqrt(T) Z),   Z ~ N(0, 1)
///
/// The estimator is the discounted sample mean of the payoff. The standard error is
/// s / sqrt(N) with s the sample standard deviation of the discounted payoffs, and the
/// 95% confidence interval is price +/- 1.96 SE (central limit theorem).
///
/// Works for any Payoff. Pricing is const and deterministic: each call re-seeds a fresh
/// generator with the stored seed, so repeated calls (for instance with bumped market
/// data) reuse the same random numbers.
class MonteCarloPricer final : public Pricer {
public:
    static constexpr std::uint64_t kDefaultSeed = 42;

    /// @param numPaths number of simulated paths, at least 2 (needed for a variance)
    /// @param seed     seed of the random number generator
    /// @throws std::invalid_argument if numPaths < 2.
    explicit MonteCarloPricer(std::size_t numPaths, std::uint64_t seed = kDefaultSeed);

    /// @throws std::invalid_argument for American exercise.
    [[nodiscard]] PricingResult price(const Option& option,
                                      const MarketData& market) const override;

    [[nodiscard]] std::size_t numPaths() const noexcept { return numPaths_; }
    [[nodiscard]] std::uint64_t seed() const noexcept { return seed_; }

    /// Two-sided 97.5% quantile of N(0,1), used for the 95% confidence interval.
    static constexpr double kZ95 = 1.959963984540054;

private:
    std::size_t numPaths_;
    std::uint64_t seed_;
};

}  // namespace pricer
