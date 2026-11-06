#pragma once

#include "pricer/greeks.hpp"
#include "pricer/pricer.hpp"

namespace pricer {

/// Closed-form Black-Scholes-Merton pricer.
///
/// Supports European options on CallPayoff, PutPayoff, DigitalCallPayoff and
/// DigitalPutPayoff. The standard error is zero (the price is exact within the model).
class AnalyticPricer final : public Pricer {
public:
    /// @throws std::invalid_argument for American exercise or an unsupported payoff.
    [[nodiscard]] PricingResult price(const Option& option,
                                      const MarketData& market) const override;

    /// Closed-form Greeks of a European call or put.
    /// @throws std::invalid_argument for American exercise or a non-vanilla payoff.
    [[nodiscard]] Greeks greeks(const Option& option, const MarketData& market) const;
};

}  // namespace pricer
