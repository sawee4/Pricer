#pragma once

namespace pricer {

enum class OptionType { Call, Put };

/// Closed-form Black-Scholes-Merton formulas with a continuous dividend yield q.
///
/// Under the risk-neutral measure S_T = F exp(-v^2/2 + v Z), Z ~ N(0,1), with
/// forward F = S exp((r - q) T) and total volatility v = sigma sqrt(T). Then
///
///     call = DF [F N(d1) - K N(d2)],   put = DF [K N(-d2) - F N(-d1)],
///     d1 = [ln(F/K) + v^2/2] / v,       d2 = d1 - v,      DF = exp(-r T).
///
/// When v -> 0 (sigma -> 0 or T -> 0) d1, d2 -> +/-inf and the formula is 0/0;
/// the limit is the discounted intrinsic value of the forward,
/// DF max(F - K, 0) for a call, which is returned explicitly below a threshold.
namespace black_scholes {

/// Below this total volatility sigma*sqrt(T) the deterministic limit is used.
inline constexpr double kMinTotalVolatility = 1e-12;

/// Inputs of a single closed-form evaluation (kept as a plain aggregate so the
/// formulas can be called directly, independently of Option/MarketData).
struct Inputs {
    double spot;        ///< S > 0
    double strike;      ///< K > 0
    double rate;        ///< r
    double dividend;    ///< q
    double volatility;  ///< sigma >= 0
    double maturity;    ///< T >= 0
};

/// d1 and d2. Requires sigma*sqrt(T) > 0.
struct D1D2 {
    double d1;
    double d2;
};
[[nodiscard]] D1D2 computeD1D2(const Inputs& in);

/// Forward price F = S exp((r - q) T).
[[nodiscard]] double forward(const Inputs& in);

/// Price of a European vanilla call or put.
[[nodiscard]] double vanillaPrice(OptionType type, const Inputs& in);

/// Price of a European cash-or-nothing digital paying 1 if S_T > K (call) or S_T < K (put):
/// DF N(d2) for the call, DF N(-d2) for the put.
[[nodiscard]] double digitalPrice(OptionType type, const Inputs& in);

}  // namespace black_scholes
}  // namespace pricer
