#pragma once

#include <memory>
#include <vector>

#include "pricer/black_scholes.hpp"
#include "pricer/payoff.hpp"

namespace pricer {

/// Payoff of a path-dependent contract.
///
/// The path holds the underlying at the M equally spaced monitoring dates
/// t_k = k T / M, k = 1..M (the spot S_0 at t = 0 is not included), so path.back() = S_T.
class PathPayoff {
public:
    virtual ~PathPayoff() = default;

    [[nodiscard]] virtual double operator()(const std::vector<double>& path) const = 0;
    [[nodiscard]] virtual std::unique_ptr<PathPayoff> clone() const = 0;

protected:
    PathPayoff() = default;
    PathPayoff(const PathPayoff&) = default;
    PathPayoff& operator=(const PathPayoff&) = default;
    PathPayoff(PathPayoff&&) = default;
    PathPayoff& operator=(PathPayoff&&) = default;
};

/// Adapter: any terminal Payoff h(S_T) seen as a path payoff (used for validation).
class TerminalPayoff final : public PathPayoff {
public:
    /// @throws std::invalid_argument if payoff is null.
    explicit TerminalPayoff(std::unique_ptr<Payoff> payoff);
    TerminalPayoff(const TerminalPayoff& other);
    TerminalPayoff& operator=(const TerminalPayoff& other);
    TerminalPayoff(TerminalPayoff&&) noexcept = default;
    TerminalPayoff& operator=(TerminalPayoff&&) noexcept = default;
    ~TerminalPayoff() override = default;

    [[nodiscard]] double operator()(const std::vector<double>& path) const override;
    [[nodiscard]] std::unique_ptr<PathPayoff> clone() const override;

private:
    std::unique_ptr<Payoff> payoff_;
};

/// Fixed-strike arithmetic Asian: (A - K)+ for a call, (K - A)+ for a put,
/// with A = (1/M) sum_k S(t_k). No closed form under Black-Scholes.
class ArithmeticAsianPayoff final : public PathPayoff {
public:
    ArithmeticAsianPayoff(OptionType type, double strike);
    [[nodiscard]] double operator()(const std::vector<double>& path) const override;
    [[nodiscard]] std::unique_ptr<PathPayoff> clone() const override;

private:
    OptionType type_;
    double strike_;
};

/// Fixed-strike geometric Asian, with G = (prod_k S(t_k))^(1/M). ln G is Gaussian under
/// Black-Scholes, so it has a closed form (closed_form::geometricAsianPrice); it is the
/// classic control variate for the arithmetic Asian (Kemna-Vorst).
class GeometricAsianPayoff final : public PathPayoff {
public:
    GeometricAsianPayoff(OptionType type, double strike);
    [[nodiscard]] double operator()(const std::vector<double>& path) const override;
    [[nodiscard]] std::unique_ptr<PathPayoff> clone() const override;

private:
    OptionType type_;
    double strike_;
};

enum class BarrierType { UpAndOut, UpAndIn, DownAndOut, DownAndIn };

/// Discretely monitored barrier option on a vanilla call or put (no rebate).
/// The barrier is hit if S(t_k) >= H (up) or S(t_k) <= H (down) at some monitoring date.
/// Knock-out: pays the vanilla payoff if never hit; knock-in: only if hit.
/// Path by path, knock-in + knock-out = vanilla.
class BarrierPayoff final : public PathPayoff {
public:
    /// @throws std::invalid_argument for a non-positive strike or barrier.
    BarrierPayoff(OptionType type, double strike, BarrierType barrierType, double barrier);
    [[nodiscard]] double operator()(const std::vector<double>& path) const override;
    [[nodiscard]] std::unique_ptr<PathPayoff> clone() const override;

private:
    OptionType type_;
    double strike_;
    BarrierType barrierType_;
    double barrier_;
};

}  // namespace pricer
