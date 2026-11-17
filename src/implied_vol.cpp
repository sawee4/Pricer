#include "pricer/implied_vol.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pricer {

ImpliedVolResult impliedVolatility(OptionType type, double price,
                                   const black_scholes::Inputs& in,
                                   const ImpliedVolSettings& settings) {
    if (!(in.maturity > 0.0) || !(in.spot > 0.0) || !(in.strike > 0.0) || !std::isfinite(price)) {
        throw std::invalid_argument("impliedVolatility: need T > 0, S > 0, K > 0, finite price");
    }

    // No-arbitrage bounds: the price must lie strictly between the sigma -> 0 limit
    // (discounted forward intrinsic) and the sigma -> infinity limit.
    const double discountedSpot = in.spot * std::exp(-in.dividend * in.maturity);
    const double discountedStrike = in.strike * std::exp(-in.rate * in.maturity);
    const double forwardIntrinsic = type == OptionType::Call ? discountedSpot - discountedStrike
                                                             : discountedStrike - discountedSpot;
    const double lower = std::max(forwardIntrinsic, 0.0);
    const double upper = type == OptionType::Call ? discountedSpot : discountedStrike;
    if (!(price > lower && price < upper)) {
        throw std::domain_error("impliedVolatility: price outside the no-arbitrage bounds");
    }

    black_scholes::Inputs trial = in;
    auto objective = [&](double sigma) {
        trial.volatility = sigma;
        return black_scholes::vanillaPrice(type, trial) - price;
    };

    double lo = settings.lowerBound;
    double hi = settings.upperBound;
    if (objective(hi) < 0.0) {
        throw std::domain_error("impliedVolatility: implied volatility above the upper bound");
    }
    if (objective(lo) > 0.0) {
        throw std::domain_error("impliedVolatility: implied volatility below the lower bound");
    }

    // Manaster-Koehler starting point: the inflection point of BS(sigma), where Newton
    // started on an increasing convex-then-concave function cannot overshoot badly.
    const double logMoneyness = std::log(black_scholes::forward(in) / in.strike);
    double sigma = std::sqrt(2.0 * std::abs(logMoneyness) / in.maturity);
    if (!(sigma > lo && sigma < hi) || sigma < 1e-3) sigma = 0.2;

    const double tolerance = settings.priceTolerance * std::max(1.0, price);
    ImpliedVolResult result;
    for (std::size_t iter = 1; iter <= settings.maxIterations; ++iter) {
        result.iterations = iter;
        const double f = objective(sigma);
        if (std::abs(f) < tolerance) {
            result.volatility = sigma;
            return result;
        }
        // Shrink the bracket: BS is increasing in sigma.
        (f > 0.0 ? hi : lo) = sigma;
        if (hi - lo < settings.volTolerance) {
            result.volatility = 0.5 * (lo + hi);
            return result;
        }

        trial.volatility = sigma;
        const double vega = black_scholes::vanillaGreeks(type, trial).vega;
        const double newton = vega > 0.0 ? sigma - f / vega : lo - 1.0;  // force bisection
        if (newton > lo && newton < hi) {
            sigma = newton;
        } else {
            sigma = 0.5 * (lo + hi);
            ++result.bisectionSteps;
        }
    }
    throw std::domain_error("impliedVolatility: did not converge");
}

}  // namespace pricer
