#pragma once

namespace pricer {

/// Market inputs of the Black-Scholes model under the risk-neutral measure:
///
///     dS_t = (r - q) S_t dt + sigma S_t dW_t
///
/// All rates are continuously compounded and annualised; volatility is annualised.
/// The object is immutable and validated at construction, so any MarketData
/// instance that exists is guaranteed to be usable by a pricer.
class MarketData {
public:
    /// @param spot       S_0 > 0
    /// @param rate       risk-free rate r (may be negative)
    /// @param dividend   continuous dividend yield q (may be negative, e.g. FX)
    /// @param volatility sigma >= 0
    /// @throws std::invalid_argument on non-finite input, S_0 <= 0 or sigma < 0.
    MarketData(double spot, double rate, double dividend, double volatility);

    [[nodiscard]] double spot() const noexcept { return spot_; }
    [[nodiscard]] double rate() const noexcept { return rate_; }
    [[nodiscard]] double dividend() const noexcept { return dividend_; }
    [[nodiscard]] double volatility() const noexcept { return volatility_; }

    /// Returns a copy with a different spot (used for bump-and-revalue Greeks).
    [[nodiscard]] MarketData withSpot(double spot) const;
    /// Returns a copy with a different volatility.
    [[nodiscard]] MarketData withVolatility(double volatility) const;
    /// Returns a copy with a different risk-free rate.
    [[nodiscard]] MarketData withRate(double rate) const;

private:
    double spot_;
    double rate_;
    double dividend_;
    double volatility_;
};

}  // namespace pricer
