#pragma once

#include <cstddef>

#include "pricer/black_scholes.hpp"

namespace pricer {

/// Result of an implied volatility search.
struct ImpliedVolResult {
    double volatility = 0.0;
    std::size_t iterations = 0;
    std::size_t bisectionSteps = 0;  ///< iterations where the Newton step was rejected
};

/// Settings of the implied volatility solver.
struct ImpliedVolSettings {
    double priceTolerance = 1e-12;   ///< stop when |BS(sigma) - price| < tol * max(1, price)
    double volTolerance = 1e-14;     ///< stop when the bracket is narrower than this
    double lowerBound = 1e-8;        ///< search interval for sigma
    double upperBound = 10.0;
    std::size_t maxIterations = 200;
};

/// Black-Scholes implied volatility of a European call or put: the sigma such that
/// black_scholes::vanillaPrice(type, {S, K, r, q, sigma, T}) = price.
///
/// Safeguarded Newton-Raphson: Newton steps sigma -= (BS(sigma) - price) / vega, inside a
/// bracket [lo, hi] that always contains the root (BS is increasing in sigma). Whenever the
/// Newton step leaves the bracket or vega is too small, a bisection step is taken instead.
/// This keeps the quadratic convergence of Newton near the root and the guaranteed
/// convergence of bisection far from it (deep ITM/OTM options, where vega is tiny).
///
/// The starting point is sigma0 = sqrt(2 |ln(F/K)| / T) (Manaster-Koehler), or 0.2 at the money.
///
/// @param in     market and contract inputs; in.volatility is ignored
/// @throws std::invalid_argument for T <= 0 or invalid inputs,
/// @throws std::domain_error if the price violates the no-arbitrage bounds
///         (no volatility can reproduce it) or the solver does not converge.
[[nodiscard]] ImpliedVolResult impliedVolatility(OptionType type, double price,
                                                 const black_scholes::Inputs& in,
                                                 const ImpliedVolSettings& settings = {});

}  // namespace pricer
