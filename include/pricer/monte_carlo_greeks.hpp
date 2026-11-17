#pragma once

#include <cstddef>
#include <cstdint>

#include "pricer/greeks.hpp"
#include "pricer/market_data.hpp"
#include "pricer/monte_carlo_pricer.hpp"
#include "pricer/option.hpp"

namespace pricer::mc_greeks {

/// A Monte Carlo estimate with its standard error.
struct Estimate {
    double value = 0.0;
    double standardError = 0.0;
};

/// Settings of the bump-and-revalue (finite difference) Greeks.
struct BumpSettings {
    double relativeSpotBump = 0.01;   ///< h = relativeSpotBump * S for delta and gamma
    double volBump = 0.001;           ///< absolute bump of sigma for vega
    double rateBump = 0.0001;         ///< absolute bump of r for rho
    double timeBump = 1.0 / 365.0;    ///< absolute bump of T (years) for theta
    /// true: every revaluation reuses the pricer's seed (common random numbers), so the
    /// Monte Carlo noise largely cancels in the differences. false: each revaluation uses
    /// a different seed (independent draws), shown for comparison only.
    bool commonRandomNumbers = true;
};

/// Greeks by central finite differences of Monte Carlo prices:
///   delta = [V(S+h) - V(S-h)] / 2h,   gamma = [V(S+h) - 2V(S) + V(S-h)] / h^2,
///   vega, rho and theta = -dV/dT by central differences in sigma, r and T.
/// Works for any payoff the pricer supports. With common random numbers the variance of
/// delta/vega/rho/theta is O(1/N) independently of h for Lipschitz payoffs; gamma of a
/// kinked payoff still has variance O(1/(h N)).
/// @throws std::invalid_argument if a bump is not positive, if T <= timeBump or if
///         sigma <= volBump (the down-bumped inputs would be invalid).
[[nodiscard]] Greeks bumpAndRevalue(const MonteCarloPricer& pricer, const Option& option,
                                    const MarketData& market, const BumpSettings& settings = {});

/// Pathwise delta: E[e^{-rT} h'(S_T) dS_T/dS_0] with dS_T/dS_0 = S_T / S_0.
/// Unbiased with no bump size, but requires a Lipschitz payoff: only European calls and puts.
/// @throws std::invalid_argument for other payoffs (e.g. digitals, whose h' is 0 a.e.),
///         American exercise, numPaths < 2 or sigma*sqrt(T) = 0.
[[nodiscard]] Estimate pathwiseDelta(const Option& option, const MarketData& market,
                                     std::size_t numPaths, std::uint64_t seed);

/// Pathwise vega: E[e^{-rT} h'(S_T) dS_T/dsigma], dS_T/dsigma = S_T (ln(S_T/S_0) - (r-q+sigma^2/2)T)/sigma.
/// Same restrictions as pathwiseDelta.
[[nodiscard]] Estimate pathwiseVega(const Option& option, const MarketData& market,
                                    std::size_t numPaths, std::uint64_t seed);

/// Likelihood ratio delta: E[e^{-rT} h(S_T) Z / (S_0 sigma sqrt(T))]. Differentiates the
/// density instead of the payoff, so it works for any payoff (digitals included), at the
/// price of a higher variance than the pathwise estimator.
/// @throws std::invalid_argument for American exercise, numPaths < 2 or sigma*sqrt(T) = 0.
[[nodiscard]] Estimate likelihoodRatioDelta(const Option& option, const MarketData& market,
                                            std::size_t numPaths, std::uint64_t seed);

}  // namespace pricer::mc_greeks
