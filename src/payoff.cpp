#include "pricer/payoff.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pricer {

StrikePayoff::StrikePayoff(double strike) : strike_(strike) {
    if (!std::isfinite(strike) || strike <= 0.0) {
        throw std::invalid_argument("StrikePayoff: strike must be finite and > 0");
    }
}

double CallPayoff::operator()(double spot) const { return std::max(spot - strike(), 0.0); }
std::unique_ptr<Payoff> CallPayoff::clone() const { return std::make_unique<CallPayoff>(*this); }

double PutPayoff::operator()(double spot) const { return std::max(strike() - spot, 0.0); }
std::unique_ptr<Payoff> PutPayoff::clone() const { return std::make_unique<PutPayoff>(*this); }

double DigitalCallPayoff::operator()(double spot) const { return spot > strike() ? 1.0 : 0.0; }
std::unique_ptr<Payoff> DigitalCallPayoff::clone() const {
    return std::make_unique<DigitalCallPayoff>(*this);
}

double DigitalPutPayoff::operator()(double spot) const { return spot < strike() ? 1.0 : 0.0; }
std::unique_ptr<Payoff> DigitalPutPayoff::clone() const {
    return std::make_unique<DigitalPutPayoff>(*this);
}

}  // namespace pricer
