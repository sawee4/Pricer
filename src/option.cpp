#include "pricer/option.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace pricer {

Option::Option(std::unique_ptr<Payoff> payoff, double maturity, ExerciseStyle exercise)
    : payoff_(std::move(payoff)), maturity_(maturity), exercise_(exercise) {
    if (!payoff_) {
        throw std::invalid_argument("Option: payoff must not be null");
    }
    if (!std::isfinite(maturity) || maturity < 0.0) {
        throw std::invalid_argument("Option: maturity must be finite and >= 0");
    }
}

Option::Option(const Option& other)
    : payoff_(other.payoff_->clone()), maturity_(other.maturity_), exercise_(other.exercise_) {}

Option& Option::operator=(const Option& other) {
    if (this != &other) {
        Option copy(other);       // copy-and-swap: strong exception guarantee
        *this = std::move(copy);
    }
    return *this;
}

}  // namespace pricer
