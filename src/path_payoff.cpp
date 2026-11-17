#include "pricer/path_payoff.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pricer {

namespace {

double vanilla(OptionType type, double underlying, double strike) {
    return type == OptionType::Call ? std::max(underlying - strike, 0.0)
                                    : std::max(strike - underlying, 0.0);
}

void requirePositive(double value, const char* what) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(what);
    }
}

}  // namespace

// --- TerminalPayoff --------------------------------------------------------------------------

TerminalPayoff::TerminalPayoff(std::unique_ptr<Payoff> payoff) : payoff_(std::move(payoff)) {
    if (!payoff_) {
        throw std::invalid_argument("TerminalPayoff: payoff must not be null");
    }
}

TerminalPayoff::TerminalPayoff(const TerminalPayoff& other) : payoff_(other.payoff_->clone()) {}

TerminalPayoff& TerminalPayoff::operator=(const TerminalPayoff& other) {
    if (this != &other) {
        payoff_ = other.payoff_->clone();
    }
    return *this;
}

double TerminalPayoff::operator()(const std::vector<double>& path) const {
    return (*payoff_)(path.back());
}

std::unique_ptr<PathPayoff> TerminalPayoff::clone() const {
    return std::make_unique<TerminalPayoff>(*this);
}

// --- Asian -----------------------------------------------------------------------------------

ArithmeticAsianPayoff::ArithmeticAsianPayoff(OptionType type, double strike)
    : type_(type), strike_(strike) {
    requirePositive(strike, "ArithmeticAsianPayoff: strike must be > 0");
}

double ArithmeticAsianPayoff::operator()(const std::vector<double>& path) const {
    double sum = 0.0;
    for (double s : path) sum += s;
    return vanilla(type_, sum / static_cast<double>(path.size()), strike_);
}

std::unique_ptr<PathPayoff> ArithmeticAsianPayoff::clone() const {
    return std::make_unique<ArithmeticAsianPayoff>(*this);
}

GeometricAsianPayoff::GeometricAsianPayoff(OptionType type, double strike)
    : type_(type), strike_(strike) {
    requirePositive(strike, "GeometricAsianPayoff: strike must be > 0");
}

double GeometricAsianPayoff::operator()(const std::vector<double>& path) const {
    // Average of logs rather than a product of prices: the product of 250 prices around 100
    // would overflow a double (100^250 = 1e500).
    double sumLog = 0.0;
    for (double s : path) sumLog += std::log(s);
    return vanilla(type_, std::exp(sumLog / static_cast<double>(path.size())), strike_);
}

std::unique_ptr<PathPayoff> GeometricAsianPayoff::clone() const {
    return std::make_unique<GeometricAsianPayoff>(*this);
}

// --- Barrier ---------------------------------------------------------------------------------

BarrierPayoff::BarrierPayoff(OptionType type, double strike, BarrierType barrierType,
                             double barrier)
    : type_(type), strike_(strike), barrierType_(barrierType), barrier_(barrier) {
    requirePositive(strike, "BarrierPayoff: strike must be > 0");
    requirePositive(barrier, "BarrierPayoff: barrier must be > 0");
}

double BarrierPayoff::operator()(const std::vector<double>& path) const {
    const bool up = barrierType_ == BarrierType::UpAndOut || barrierType_ == BarrierType::UpAndIn;
    const bool knockIn =
        barrierType_ == BarrierType::UpAndIn || barrierType_ == BarrierType::DownAndIn;
    const bool hit = up ? std::any_of(path.begin(), path.end(),
                                      [this](double s) { return s >= barrier_; })
                        : std::any_of(path.begin(), path.end(),
                                      [this](double s) { return s <= barrier_; });
    const bool alive = knockIn ? hit : !hit;
    return alive ? vanilla(type_, path.back(), strike_) : 0.0;
}

std::unique_ptr<PathPayoff> BarrierPayoff::clone() const {
    return std::make_unique<BarrierPayoff>(*this);
}

}  // namespace pricer
