#pragma once

#include <cmath>

namespace pricer::math {

/// 1 / sqrt(2*pi)
inline constexpr double kInvSqrt2Pi = 0.39894228040143267794;
/// 1 / sqrt(2)
inline constexpr double kInvSqrt2 = 0.70710678118654752440;

/// Standard normal density n(x) = exp(-x^2/2) / sqrt(2*pi).
[[nodiscard]] inline double normPdf(double x) noexcept {
    return kInvSqrt2Pi * std::exp(-0.5 * x * x);
}

/// Standard normal CDF N(x) = P(Z <= x).
///
/// Written as erfc(-x/sqrt 2)/2 rather than (1 + erf(x/sqrt 2))/2: in the left
/// tail the latter computes 1 + (-1 + tiny), which cancels to 0 well before the
/// true value underflows (e.g. N(-10) ~ 7.6e-24 would come out as exactly 0).
/// erfc keeps full relative precision in the tail and never overflows.
[[nodiscard]] inline double normCdf(double x) noexcept {
    return 0.5 * std::erfc(-x * kInvSqrt2);
}

}  // namespace pricer::math
