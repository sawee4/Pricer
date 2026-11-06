#pragma once

#include <memory>

#include "pricer/payoff.hpp"

namespace pricer {

enum class ExerciseStyle {
    European,  ///< exercise at maturity only
    American   ///< exercise at any time in [0, T]
};

/// A contract = payoff + maturity + exercise style.
///
/// Owns its payoff (RAII through std::unique_ptr) and has value semantics:
/// copying an Option deep-copies the payoff via Payoff::clone().
class Option {
public:
    /// @param payoff   terminal payoff, must not be null
    /// @param maturity time to expiry T in years, T >= 0
    /// @throws std::invalid_argument on null payoff or invalid maturity.
    Option(std::unique_ptr<Payoff> payoff, double maturity,
           ExerciseStyle exercise = ExerciseStyle::European);

    Option(const Option& other);
    Option& operator=(const Option& other);
    Option(Option&&) noexcept = default;
    Option& operator=(Option&&) noexcept = default;
    ~Option() = default;

    [[nodiscard]] const Payoff& payoff() const noexcept { return *payoff_; }
    [[nodiscard]] double maturity() const noexcept { return maturity_; }
    [[nodiscard]] ExerciseStyle exercise() const noexcept { return exercise_; }

    /// Convenience: payoff evaluated at a given underlying price.
    [[nodiscard]] double payoffAt(double spot) const { return (*payoff_)(spot); }

private:
    std::unique_ptr<Payoff> payoff_;
    double maturity_;
    ExerciseStyle exercise_;
};

}  // namespace pricer
