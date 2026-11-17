#pragma once

#include <cstddef>

#include "pricer/pricer.hpp"

namespace pricer {

/// Cox-Ross-Rubinstein binomial tree for European and American options.
///
/// With n steps of length dt = T / n, the underlying moves up by u = exp(sigma sqrt(dt))
/// or down by d = 1 / u. The risk-neutral up probability
///
///     p = (exp((r - q) dt) - d) / (u - d)
///
/// makes the discounted price a martingale. Prices are obtained by backward induction:
/// V = exp(-r dt) [p V_up + (1 - p) V_down] (continuation), and for American exercise
/// V = max(continuation, payoff(S)) at every node (dynamic programming).
///
/// European prices converge to Black-Scholes in O(1/n), with an oscillation that depends
/// on where the strike falls between terminal nodes. Cost is O(n^2) time, O(n) memory.
class BinomialPricer final : public Pricer {
public:
    /// @throws std::invalid_argument if numSteps == 0.
    explicit BinomialPricer(std::size_t numSteps);

    /// @throws std::invalid_argument if sigma = 0 with T > 0 (the tree degenerates),
    /// @throws std::domain_error if p is outside [0, 1] (no-arbitrage violated: increase
    ///         the number of steps).
    [[nodiscard]] PricingResult price(const Option& option,
                                      const MarketData& market) const override;

    [[nodiscard]] std::size_t numSteps() const noexcept { return numSteps_; }

private:
    std::size_t numSteps_;
};

}  // namespace pricer
