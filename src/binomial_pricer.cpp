#include "pricer/binomial_pricer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace pricer {

BinomialPricer::BinomialPricer(std::size_t numSteps) : numSteps_(numSteps) {
    if (numSteps == 0) {
        throw std::invalid_argument("BinomialPricer: need at least one step");
    }
}

PricingResult BinomialPricer::price(const Option& option, const MarketData& market) const {
    const Payoff& payoff = option.payoff();
    const double t = option.maturity();
    const double spot = market.spot();

    if (t == 0.0) {
        const double intrinsic = payoff(spot);
        return PricingResult{intrinsic, 0.0, intrinsic, intrinsic};
    }
    if (!(market.volatility() > 0.0)) {
        throw std::invalid_argument("BinomialPricer: volatility must be > 0");
    }

    const std::size_t n = numSteps_;
    const double dt = t / static_cast<double>(n);
    const double logU = market.volatility() * std::sqrt(dt);  // ln u; ln d = -ln u
    const double u = std::exp(logU);
    const double d = 1.0 / u;
    const double growth = std::exp((market.rate() - market.dividend()) * dt);
    const double p = (growth - d) / (u - d);
    if (!(p >= 0.0 && p <= 1.0)) {
        throw std::domain_error(
            "BinomialPricer: risk-neutral probability outside [0, 1] (arbitrage in the tree); "
            "increase the number of steps");
    }
    const double discount = std::exp(-market.rate() * dt);
    const double discP = discount * p;
    const double discQ = discount * (1.0 - p);
    const bool american = option.exercise() == ExerciseStyle::American;

    // Node j at step i has S = S_0 u^j d^(i-j) = S_0 exp((2j - i) ln u). Computing it directly
    // (rather than multiplying by u repeatedly) avoids accumulating rounding errors.
    auto nodeSpot = [&](std::size_t i, std::size_t j) {
        return spot * std::exp((2.0 * static_cast<double>(j) - static_cast<double>(i)) * logU);
    };

    // A single vector updated in place: values[j] at step i only needs values[j] and
    // values[j + 1] from step i + 1, and j increases, so nothing is overwritten too early.
    std::vector<double> values(n + 1);
    for (std::size_t j = 0; j <= n; ++j) {
        values[j] = payoff(nodeSpot(n, j));
    }
    for (std::size_t i = n; i-- > 0;) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double continuation = discP * values[j + 1] + discQ * values[j];
            values[j] = american ? std::max(continuation, payoff(nodeSpot(i, j))) : continuation;
        }
    }

    const double value = values[0];
    return PricingResult{value, 0.0, value, value};
}

}  // namespace pricer
