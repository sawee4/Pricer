#pragma once

#include "pricer/market_data.hpp"
#include "pricer/option.hpp"

namespace pricer {

/// Result of a pricing call. Deterministic pricers (analytic, binomial)
/// report a zero standard error; Monte Carlo pricers fill all fields.
struct PricingResult {
    double price = 0.0;
    double standardError = 0.0;
    double confidenceLow = 0.0;   ///< lower bound of the 95% confidence interval
    double confidenceHigh = 0.0;  ///< upper bound of the 95% confidence interval
};

/// Common interface of all pricing engines.
///
/// A Pricer holds only its numerical settings (number of paths, steps, seed...);
/// the contract and the market are passed at call time, so the same engine
/// can price many options and pricing is a const operation.
class Pricer {
public:
    virtual ~Pricer() = default;

    [[nodiscard]] virtual PricingResult price(const Option& option,
                                              const MarketData& market) const = 0;

protected:
    Pricer() = default;
    Pricer(const Pricer&) = default;
    Pricer& operator=(const Pricer&) = default;
    Pricer(Pricer&&) = default;
    Pricer& operator=(Pricer&&) = default;
};

}  // namespace pricer
