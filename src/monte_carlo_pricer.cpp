#include "pricer/monte_carlo_pricer.hpp"

#include <cmath>
#include <stdexcept>

#include "pricer/random.hpp"
#include "pricer/statistics.hpp"

namespace pricer {

namespace {

/// Exact GBM terminal value S_T = S_0 exp(drift + diffusion * z), with constants hoisted.
struct TerminalSampler {
    double spot;
    double drift;      ///< (r - q - sigma^2/2) T
    double diffusion;  ///< sigma sqrt(T)

    [[nodiscard]] double operator()(double z) const { return spot * std::exp(drift + diffusion * z); }
};

PricingResult makeResult(double price, double variance, std::size_t samples) {
    const double se = std::sqrt(variance / static_cast<double>(samples));
    return PricingResult{price, se, price - MonteCarloPricer::kZ95 * se,
                         price + MonteCarloPricer::kZ95 * se};
}

}  // namespace

MonteCarloPricer::MonteCarloPricer(std::size_t numPaths, std::uint64_t seed,
                                   VarianceReduction reduction)
    : numPaths_(numPaths), seed_(seed), reduction_(reduction) {
    const std::size_t minimum = reduction == VarianceReduction::Antithetic ? 4 : 2;
    if (numPaths < minimum) {
        throw std::invalid_argument("MonteCarloPricer: too few paths for a variance estimate");
    }
}

PricingResult MonteCarloPricer::price(const Option& option, const MarketData& market) const {
    if (option.exercise() != ExerciseStyle::European) {
        throw std::invalid_argument(
            "MonteCarloPricer: only European exercise supported (use BinomialPricer)");
    }

    const double t = option.maturity();
    const double vol = market.volatility();
    const TerminalSampler terminal{market.spot(),
                                   (market.rate() - market.dividend() - 0.5 * vol * vol) * t,
                                   vol * std::sqrt(t)};
    const double discount = std::exp(-market.rate() * t);
    const Payoff& payoff = option.payoff();
    NormalGenerator rng(seed_);

    switch (reduction_) {
        case VarianceReduction::None: {
            RunningStats stats;
            for (std::size_t i = 0; i < numPaths_; ++i) {
                stats.add(discount * payoff(terminal(rng.nextNormal())));
            }
            return makeResult(stats.mean(), stats.variance(), numPaths_);
        }

        case VarianceReduction::Antithetic: {
            // The two halves of a pair are dependent, so the i.i.d. samples are the pair
            // averages: the variance and SE must be computed on them, not on 2N values.
            const std::size_t pairs = numPaths_ / 2;
            RunningStats stats;
            for (std::size_t i = 0; i < pairs; ++i) {
                const double z = rng.nextNormal();
                stats.add(0.5 * discount * (payoff(terminal(z)) + payoff(terminal(-z))));
            }
            return makeResult(stats.mean(), stats.variance(), pairs);
        }

        case VarianceReduction::ControlVariate: {
            // Control X = e^{-rT} S_T, E[X] = S_0 e^{-qT} (the discounted forward).
            const double expectedControl = market.spot() * std::exp(-market.dividend() * t);
            RunningCovariance stats;  // X = control, Y = discounted payoff
            for (std::size_t i = 0; i < numPaths_; ++i) {
                const double sT = terminal(rng.nextNormal());
                stats.add(discount * sT, discount * payoff(sT));
            }
            const double varX = stats.varianceX();
            // sigma = 0 or T = 0: the control is constant and carries no information.
            const double beta = varX > 0.0 ? stats.covariance() / varX : 0.0;
            const double price = stats.meanY() - beta * (stats.meanX() - expectedControl);
            // Var(Y - beta X) = Var(Y) - 2 beta Cov + beta^2 Var(X) = Var(Y) (1 - rho^2) at beta*.
            double variance = stats.varianceY() - 2.0 * beta * stats.covariance() +
                              beta * beta * varX;
            variance = variance > 0.0 ? variance : 0.0;  // guard against rounding below zero
            return makeResult(price, variance, numPaths_);
        }
    }
    throw std::logic_error("MonteCarloPricer: unknown variance reduction");
}

}  // namespace pricer
