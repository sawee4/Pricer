#pragma once

#include <cstddef>

#include "pricer/black_scholes.hpp"

/// Closed-form prices of exotic options under Black-Scholes, used to validate the
/// path-dependent Monte Carlo pricer.
namespace pricer::closed_form {

/// Discretely monitored fixed-strike geometric Asian option, fixings at t_k = k T / M, k = 1..M.
///
/// ln G = ln S_0 + (1/M) sum_k ln(S(t_k)/S_0) is Gaussian with
///   mean     mu   = ln S_0 + (r - q - sigma^2/2) T (M + 1) / (2M)
///   variance s^2  = sigma^2 T (M + 1)(2M + 1) / (6 M^2)
/// so the price is a Black-type formula: e^{-rT} [e^{mu + s^2/2} N(d1) - K N(d2)] for a call,
/// d1 = (mu - ln K + s^2) / s, d2 = d1 - s. For M = 1 it reduces to Black-Scholes.
/// @throws std::invalid_argument if numFixings == 0, T <= 0 or sigma <= 0.
[[nodiscard]] double geometricAsianPrice(OptionType type, const black_scholes::Inputs& in,
                                         std::size_t numFixings);

/// Continuously monitored down-and-in call with barrier H <= K and H < S (Hull, ch. 26):
///   c_di = S e^{-qT} (H/S)^{2 lambda} N(y) - K e^{-rT} (H/S)^{2 lambda - 2} N(y - sigma sqrt T),
///   lambda = (r - q + sigma^2/2) / sigma^2,
///   y = ln(H^2 / (S K)) / (sigma sqrt T) + lambda sigma sqrt T.
/// @throws std::invalid_argument unless 0 < H <= K, H < S, T > 0 and sigma > 0.
[[nodiscard]] double downAndInCallPrice(const black_scholes::Inputs& in, double barrier);

/// Continuously monitored down-and-out call, H <= K: vanilla call - down-and-in call.
[[nodiscard]] double downAndOutCallPrice(const black_scholes::Inputs& in, double barrier);

/// Broadie-Glasserman-Kou continuity correction: a barrier monitored at intervals dt behaves
/// like a continuously monitored barrier shifted away from the spot,
///   H' = H exp(+beta sigma sqrt(dt)) for an up barrier, H exp(-beta sigma sqrt(dt)) for a down
///   barrier, with beta = -zeta(1/2) / sqrt(2 pi) ~ 0.5826.
inline constexpr double kBgkBeta = 0.5825971579390106;
[[nodiscard]] double continuityCorrectedBarrier(double barrier, double volatility, double dt,
                                                bool upBarrier);

}  // namespace pricer::closed_form
