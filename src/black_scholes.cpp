#include "pricer/black_scholes.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "pricer/normal.hpp"

namespace pricer::black_scholes {

namespace {

void validate(const Inputs& in) {
    if (!(in.spot > 0.0) || !(in.strike > 0.0) || !(in.volatility >= 0.0) ||
        !(in.maturity >= 0.0) || !std::isfinite(in.rate) || !std::isfinite(in.dividend) ||
        !std::isfinite(in.spot) || !std::isfinite(in.strike) || !std::isfinite(in.volatility) ||
        !std::isfinite(in.maturity)) {
        throw std::invalid_argument("black_scholes: invalid inputs");
    }
}

double totalVolatility(const Inputs& in) { return in.volatility * std::sqrt(in.maturity); }

}  // namespace

double forward(const Inputs& in) {
    return in.spot * std::exp((in.rate - in.dividend) * in.maturity);
}

D1D2 computeD1D2(const Inputs& in) {
    const double v = totalVolatility(in);
    if (!(v > 0.0)) {
        throw std::domain_error("black_scholes: d1/d2 undefined for zero total volatility");
    }
    // ln(F/K) = ln(S/K) + (r - q) T, written via the forward to keep the formula compact.
    const double d1 = (std::log(forward(in) / in.strike) + 0.5 * v * v) / v;
    return {d1, d1 - v};
}

double vanillaPrice(OptionType type, const Inputs& in) {
    validate(in);
    const double df = std::exp(-in.rate * in.maturity);
    const double fwd = forward(in);
    const double k = in.strike;

    if (totalVolatility(in) < kMinTotalVolatility) {
        // Deterministic limit: S_T = F almost surely.
        return type == OptionType::Call ? df * std::max(fwd - k, 0.0)
                                        : df * std::max(k - fwd, 0.0);
    }

    const auto [d1, d2] = computeD1D2(in);
    if (type == OptionType::Call) {
        return df * (fwd * math::normCdf(d1) - k * math::normCdf(d2));
    }
    return df * (k * math::normCdf(-d2) - fwd * math::normCdf(-d1));
}

double digitalPrice(OptionType type, const Inputs& in) {
    validate(in);
    const double df = std::exp(-in.rate * in.maturity);
    const double fwd = forward(in);

    if (totalVolatility(in) < kMinTotalVolatility) {
        // Same strict-inequality convention as DigitalCallPayoff / DigitalPutPayoff.
        const bool inTheMoney = type == OptionType::Call ? fwd > in.strike : fwd < in.strike;
        return inTheMoney ? df : 0.0;
    }

    const double d2 = computeD1D2(in).d2;
    return df * math::normCdf(type == OptionType::Call ? d2 : -d2);
}

}  // namespace pricer::black_scholes
