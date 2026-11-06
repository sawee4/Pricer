#pragma once

#include <memory>

namespace pricer {

/// Terminal payoff h(S_T) of a contract.
///
/// Payoffs are polymorphic and immutable. They are owned through
/// std::unique_ptr and copied with clone() ("virtual copy constructor"),
/// so that an Option can hold any payoff by value semantics.
class Payoff {
public:
    virtual ~Payoff() = default;

    /// Payoff amount for a given underlying price at expiry.
    [[nodiscard]] virtual double operator()(double spot) const = 0;

    /// Deep copy preserving the dynamic type.
    [[nodiscard]] virtual std::unique_ptr<Payoff> clone() const = 0;

protected:
    Payoff() = default;
    Payoff(const Payoff&) = default;
    Payoff& operator=(const Payoff&) = default;
    Payoff(Payoff&&) = default;
    Payoff& operator=(Payoff&&) = default;
};

/// Base class for payoffs parameterised by a single strike K > 0.
class StrikePayoff : public Payoff {
public:
    /// @throws std::invalid_argument if strike is not strictly positive and finite.
    explicit StrikePayoff(double strike);

    [[nodiscard]] double strike() const noexcept { return strike_; }

private:
    double strike_;
};

/// Vanilla call: max(S - K, 0).
class CallPayoff final : public StrikePayoff {
public:
    using StrikePayoff::StrikePayoff;
    [[nodiscard]] double operator()(double spot) const override;
    [[nodiscard]] std::unique_ptr<Payoff> clone() const override;
};

/// Vanilla put: max(K - S, 0).
class PutPayoff final : public StrikePayoff {
public:
    using StrikePayoff::StrikePayoff;
    [[nodiscard]] double operator()(double spot) const override;
    [[nodiscard]] std::unique_ptr<Payoff> clone() const override;
};

/// Cash-or-nothing digital call: 1 if S > K, 0 otherwise.
class DigitalCallPayoff final : public StrikePayoff {
public:
    using StrikePayoff::StrikePayoff;
    [[nodiscard]] double operator()(double spot) const override;
    [[nodiscard]] std::unique_ptr<Payoff> clone() const override;
};

/// Cash-or-nothing digital put: 1 if S < K, 0 otherwise.
class DigitalPutPayoff final : public StrikePayoff {
public:
    using StrikePayoff::StrikePayoff;
    [[nodiscard]] double operator()(double spot) const override;
    [[nodiscard]] std::unique_ptr<Payoff> clone() const override;
};

}  // namespace pricer
