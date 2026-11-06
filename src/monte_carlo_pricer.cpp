#include "pricer/monte_carlo_pricer.hpp"

#include <cmath>
#include <stdexcept>

#include "pricer/random.hpp"
#include "pricer/statistics.hpp"

namespace pricer {

MonteCarloPricer::MonteCarloPricer(std::size_t numPaths, std::uint64_t seed)
    : numPaths_(numPaths), seed_(seed) {
    if (numPaths < 2) {
        throw std::invalid_argument("MonteCarloPricer: need at least 2 paths");
    }
}

PricingResult MonteCarloPricer::price(const Option& option, const MarketData& market) const {
    if (option.exercise() != ExerciseStyle::European) {
        throw std::invalid_argument(
            "MonteCarloPricer: only European exercise supported (use BinomialPricer)");
    }

    const double t = option.maturity();
    const double vol = market.volatility();
    // Hoist everything that does not depend on the draw out of the loop.
    const double drift = (market.rate() - market.dividend() - 0.5 * vol * vol) * t;
    const double diffusion = vol * std::sqrt(t);
    const double discount = std::exp(-market.rate() * t);
    const double spot = market.spot();
    const Payoff& payoff = option.payoff();

    NormalGenerator rng(seed_);
    RunningStats stats;
    for (std::size_t i = 0; i < numPaths_; ++i) {
        const double z = rng.nextNormal();
        const double terminal = spot * std::exp(drift + diffusion * z);
        stats.add(discount * payoff(terminal));
    }

    const double price = stats.mean();
    const double standardError = std::sqrt(stats.variance() / static_cast<double>(numPaths_));
    return PricingResult{price, standardError, price - kZ95 * standardError,
                         price + kZ95 * standardError};
}

}  // namespace pricer
