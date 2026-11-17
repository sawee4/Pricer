#include "pricer/path_dependent_pricer.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "pricer/monte_carlo_pricer.hpp"
#include "pricer/random.hpp"
#include "pricer/statistics.hpp"

namespace pricer {

namespace {

/// Generates GBM paths at equally spaced dates into a reusable buffer (no allocation per path).
class PathGenerator {
public:
    PathGenerator(const MarketData& market, double maturity, std::size_t steps, std::uint64_t seed)
        : rng_(seed), path_(steps), spot_(market.spot()) {
        const double dt = maturity / static_cast<double>(steps);
        const double vol = market.volatility();
        drift_ = (market.rate() - market.dividend() - 0.5 * vol * vol) * dt;
        diffusion_ = vol * std::sqrt(dt);
    }

    const std::vector<double>& next() {
        // Accumulate in log space: one exp per step, and no drift of rounding errors from
        // repeated multiplications.
        double logS = std::log(spot_);
        for (double& s : path_) {
            logS += drift_ + diffusion_ * rng_.nextNormal();
            s = std::exp(logS);
        }
        return path_;
    }

private:
    NormalGenerator rng_;
    std::vector<double> path_;
    double spot_;
    double drift_ = 0.0;
    double diffusion_ = 0.0;
};

PricingResult makeResult(double price, double variance, std::size_t samples) {
    const double se = std::sqrt(variance / static_cast<double>(samples));
    return PricingResult{price, se, price - MonteCarloPricer::kZ95 * se,
                         price + MonteCarloPricer::kZ95 * se};
}

void checkMaturity(double maturity) {
    if (!std::isfinite(maturity) || maturity <= 0.0) {
        throw std::invalid_argument("PathDependentPricer: maturity must be > 0");
    }
}

}  // namespace

PathDependentPricer::PathDependentPricer(std::size_t numPaths, std::size_t numSteps,
                                         std::uint64_t seed)
    : numPaths_(numPaths), numSteps_(numSteps), seed_(seed) {
    if (numPaths < 2 || numSteps == 0) {
        throw std::invalid_argument("PathDependentPricer: need numPaths >= 2 and numSteps >= 1");
    }
}

PricingResult PathDependentPricer::price(const PathPayoff& payoff, double maturity,
                                         const MarketData& market) const {
    checkMaturity(maturity);
    const double discount = std::exp(-market.rate() * maturity);
    PathGenerator generator(market, maturity, numSteps_, seed_);
    RunningStats stats;
    for (std::size_t i = 0; i < numPaths_; ++i) {
        stats.add(discount * payoff(generator.next()));
    }
    return makeResult(stats.mean(), stats.variance(), numPaths_);
}

PricingResult PathDependentPricer::priceWithControlVariate(const PathPayoff& payoff,
                                                           const PathPayoff& control,
                                                           double controlPrice, double maturity,
                                                           const MarketData& market) const {
    checkMaturity(maturity);
    const double discount = std::exp(-market.rate() * maturity);
    PathGenerator generator(market, maturity, numSteps_, seed_);
    RunningCovariance stats;  // X = control, Y = target
    for (std::size_t i = 0; i < numPaths_; ++i) {
        const std::vector<double>& path = generator.next();
        stats.add(discount * control(path), discount * payoff(path));
    }
    const double varX = stats.varianceX();
    const double beta = varX > 0.0 ? stats.covariance() / varX : 0.0;
    const double price = stats.meanY() - beta * (stats.meanX() - controlPrice);
    double variance = stats.varianceY() - 2.0 * beta * stats.covariance() + beta * beta * varX;
    variance = variance > 0.0 ? variance : 0.0;
    return makeResult(price, variance, numPaths_);
}

}  // namespace pricer
