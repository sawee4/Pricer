#pragma once

#include <cstddef>
#include <cstdint>

#include "pricer/market_data.hpp"
#include "pricer/path_payoff.hpp"
#include "pricer/pricer.hpp"

namespace pricer {

/// Monte Carlo pricer for path-dependent European payoffs under Black-Scholes.
///
/// Simulates the underlying at M equally spaced dates with the exact GBM transition
///
///     S(t_{k+1}) = S(t_k) exp((r - q - sigma^2/2) dt + sigma sqrt(dt) Z_k),  dt = T / M,
///
/// so the values at the monitoring dates have no discretisation bias. (A *continuously*
/// monitored contract priced with M dates is still biased: see closed_form.)
///
/// It does not implement the Pricer interface, whose Option only holds terminal payoffs.
class PathDependentPricer {
public:
    /// @throws std::invalid_argument if numPaths < 2 or numSteps == 0.
    PathDependentPricer(std::size_t numPaths, std::size_t numSteps, std::uint64_t seed = 42);

    /// Plain Monte Carlo estimate with standard error and 95% confidence interval.
    /// @throws std::invalid_argument if maturity <= 0.
    [[nodiscard]] PricingResult price(const PathPayoff& payoff, double maturity,
                                      const MarketData& market) const;

    /// Control variate estimate: mean(h) - beta (mean(c) - controlPrice), with c the
    /// discounted control payoff, controlPrice its exact value, and beta = Cov(h, c) / Var(c)
    /// estimated on the same paths.
    [[nodiscard]] PricingResult priceWithControlVariate(const PathPayoff& payoff,
                                                        const PathPayoff& control,
                                                        double controlPrice, double maturity,
                                                        const MarketData& market) const;

    [[nodiscard]] std::size_t numPaths() const noexcept { return numPaths_; }
    [[nodiscard]] std::size_t numSteps() const noexcept { return numSteps_; }
    [[nodiscard]] std::uint64_t seed() const noexcept { return seed_; }

private:
    std::size_t numPaths_;
    std::size_t numSteps_;
    std::uint64_t seed_;
};

}  // namespace pricer
