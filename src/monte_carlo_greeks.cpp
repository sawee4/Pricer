#include "pricer/monte_carlo_greeks.hpp"

#include <cmath>
#include <stdexcept>

#include "pricer/payoff.hpp"
#include "pricer/random.hpp"
#include "pricer/statistics.hpp"

namespace pricer::mc_greeks {

namespace {

Option withMaturity(const Option& option, double maturity) {
    return Option(option.payoff().clone(), maturity, option.exercise());
}

/// Common setup of the estimators that simulate S_T = S_0 exp(drift + diffusion * Z).
struct Simulation {
    double spot;
    double drift;
    double diffusion;
    double discount;

    Simulation(const Option& option, const MarketData& market, std::size_t numPaths) {
        if (option.exercise() != ExerciseStyle::European) {
            throw std::invalid_argument("mc_greeks: only European exercise supported");
        }
        if (numPaths < 2) {
            throw std::invalid_argument("mc_greeks: need at least 2 paths");
        }
        const double t = option.maturity();
        const double vol = market.volatility();
        spot = market.spot();
        drift = (market.rate() - market.dividend() - 0.5 * vol * vol) * t;
        diffusion = vol * std::sqrt(t);
        discount = std::exp(-market.rate() * t);
        if (!(diffusion > 0.0)) {
            throw std::invalid_argument("mc_greeks: sigma * sqrt(T) must be > 0");
        }
    }

    [[nodiscard]] double terminal(double z) const { return spot * std::exp(drift + diffusion * z); }
};

Estimate toEstimate(const RunningStats& stats) {
    return {stats.mean(), std::sqrt(stats.variance() / static_cast<double>(stats.count()))};
}

/// Returns +1 for a call, -1 for a put; throws for any other payoff.
double vanillaSign(const Payoff& payoff, double& strike) {
    if (const auto* call = dynamic_cast<const CallPayoff*>(&payoff)) {
        strike = call->strike();
        return 1.0;
    }
    if (const auto* put = dynamic_cast<const PutPayoff*>(&payoff)) {
        strike = put->strike();
        return -1.0;
    }
    throw std::invalid_argument(
        "mc_greeks: pathwise estimators need a Lipschitz payoff (vanilla call or put); "
        "use likelihoodRatioDelta for digitals");
}

}  // namespace

Greeks bumpAndRevalue(const MonteCarloPricer& pricer, const Option& option,
                      const MarketData& market, const BumpSettings& settings) {
    if (!(settings.relativeSpotBump > 0.0) || !(settings.volBump > 0.0) ||
        !(settings.rateBump > 0.0) || !(settings.timeBump > 0.0)) {
        throw std::invalid_argument("bumpAndRevalue: bumps must be > 0");
    }
    if (!(option.maturity() > settings.timeBump)) {
        throw std::invalid_argument("bumpAndRevalue: maturity must exceed the time bump");
    }
    if (!(market.volatility() > settings.volBump)) {
        throw std::invalid_argument("bumpAndRevalue: volatility must exceed the vol bump");
    }

    // Each revaluation gets its own pricer: same seed (CRN) or a distinct one (independent).
    std::uint64_t nextSeed = pricer.seed();
    auto value = [&](const Option& opt, const MarketData& mkt) {
        const std::uint64_t seed = settings.commonRandomNumbers ? pricer.seed() : nextSeed++;
        return MonteCarloPricer(pricer.numPaths(), seed, pricer.reduction()).price(opt, mkt).price;
    };

    const double s = market.spot();
    const double h = settings.relativeSpotBump * s;
    const double dv = settings.volBump;
    const double dr = settings.rateBump;
    const double dt = settings.timeBump;
    const double t = option.maturity();

    const double base = value(option, market);
    const double spotUp = value(option, market.withSpot(s + h));
    const double spotDown = value(option, market.withSpot(s - h));

    Greeks g;
    g.delta = (spotUp - spotDown) / (2.0 * h);
    g.gamma = (spotUp - 2.0 * base + spotDown) / (h * h);
    g.vega = (value(option, market.withVolatility(market.volatility() + dv)) -
              value(option, market.withVolatility(market.volatility() - dv))) / (2.0 * dv);
    g.rho = (value(option, market.withRate(market.rate() + dr)) -
             value(option, market.withRate(market.rate() - dr))) / (2.0 * dr);
    g.theta = -(value(withMaturity(option, t + dt), market) -
                value(withMaturity(option, t - dt), market)) / (2.0 * dt);
    return g;
}

Estimate pathwiseDelta(const Option& option, const MarketData& market, std::size_t numPaths,
                       std::uint64_t seed) {
    double strike = 0.0;
    const double sign = vanillaSign(option.payoff(), strike);
    const Simulation sim(option, market, numPaths);

    NormalGenerator rng(seed);
    RunningStats stats;
    for (std::size_t i = 0; i < numPaths; ++i) {
        const double sT = sim.terminal(rng.nextNormal());
        // h'(S_T) = sign * 1{sign * (S_T - K) > 0}, dS_T/dS_0 = S_T / S_0.
        const bool inTheMoney = sign * (sT - strike) > 0.0;
        stats.add(inTheMoney ? sim.discount * sign * sT / sim.spot : 0.0);
    }
    return toEstimate(stats);
}

Estimate pathwiseVega(const Option& option, const MarketData& market, std::size_t numPaths,
                      std::uint64_t seed) {
    double strike = 0.0;
    const double sign = vanillaSign(option.payoff(), strike);
    const Simulation sim(option, market, numPaths);
    const double vol = market.volatility();
    const double t = option.maturity();
    const double driftPlus = (market.rate() - market.dividend() + 0.5 * vol * vol) * t;

    NormalGenerator rng(seed);
    RunningStats stats;
    for (std::size_t i = 0; i < numPaths; ++i) {
        const double sT = sim.terminal(rng.nextNormal());
        const bool inTheMoney = sign * (sT - strike) > 0.0;
        // dS_T/dsigma = S_T (sqrt(T) Z - sigma T) = S_T (ln(S_T/S_0) - (r - q + sigma^2/2) T) / sigma
        const double dSdSigma = sT * (std::log(sT / sim.spot) - driftPlus) / vol;
        stats.add(inTheMoney ? sim.discount * sign * dSdSigma : 0.0);
    }
    return toEstimate(stats);
}

Estimate likelihoodRatioDelta(const Option& option, const MarketData& market,
                              std::size_t numPaths, std::uint64_t seed) {
    const Simulation sim(option, market, numPaths);
    const Payoff& payoff = option.payoff();
    // Score of the lognormal density with respect to S_0: Z / (S_0 sigma sqrt(T)).
    const double scale = sim.discount / (sim.spot * sim.diffusion);

    NormalGenerator rng(seed);
    RunningStats stats;
    for (std::size_t i = 0; i < numPaths; ++i) {
        const double z = rng.nextNormal();
        stats.add(scale * payoff(sim.terminal(z)) * z);
    }
    return toEstimate(stats);
}

}  // namespace pricer::mc_greeks
