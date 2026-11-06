#include "pricer/analytic_pricer.hpp"

#include <stdexcept>

#include "pricer/black_scholes.hpp"
#include "pricer/payoff.hpp"

namespace pricer {

namespace {

black_scholes::Inputs makeInputs(double strike, const Option& option, const MarketData& market) {
    return {market.spot(),     strike,           market.rate(),
            market.dividend(), market.volatility(), option.maturity()};
}

}  // namespace

PricingResult AnalyticPricer::price(const Option& option, const MarketData& market) const {
    if (option.exercise() != ExerciseStyle::European) {
        throw std::invalid_argument(
            "AnalyticPricer: only European exercise has a closed form (use BinomialPricer)");
    }

    // Dispatch on the dynamic payoff type. A closed form exists only for a known,
    // finite set of payoffs, so an explicit type check is clearer than adding a
    // pricing method to every Payoff (which would couple payoffs to one model).
    const Payoff& payoff = option.payoff();
    double value = 0.0;
    if (const auto* call = dynamic_cast<const CallPayoff*>(&payoff)) {
        value = black_scholes::vanillaPrice(OptionType::Call,
                                            makeInputs(call->strike(), option, market));
    } else if (const auto* put = dynamic_cast<const PutPayoff*>(&payoff)) {
        value = black_scholes::vanillaPrice(OptionType::Put,
                                            makeInputs(put->strike(), option, market));
    } else if (const auto* dcall = dynamic_cast<const DigitalCallPayoff*>(&payoff)) {
        value = black_scholes::digitalPrice(OptionType::Call,
                                            makeInputs(dcall->strike(), option, market));
    } else if (const auto* dput = dynamic_cast<const DigitalPutPayoff*>(&payoff)) {
        value = black_scholes::digitalPrice(OptionType::Put,
                                            makeInputs(dput->strike(), option, market));
    } else {
        throw std::invalid_argument("AnalyticPricer: unsupported payoff type");
    }

    return PricingResult{value, 0.0, value, value};
}

Greeks AnalyticPricer::greeks(const Option& option, const MarketData& market) const {
    if (option.exercise() != ExerciseStyle::European) {
        throw std::invalid_argument("AnalyticPricer::greeks: only European exercise supported");
    }
    const Payoff& payoff = option.payoff();
    if (const auto* call = dynamic_cast<const CallPayoff*>(&payoff)) {
        return black_scholes::vanillaGreeks(OptionType::Call,
                                            makeInputs(call->strike(), option, market));
    }
    if (const auto* put = dynamic_cast<const PutPayoff*>(&payoff)) {
        return black_scholes::vanillaGreeks(OptionType::Put,
                                            makeInputs(put->strike(), option, market));
    }
    throw std::invalid_argument("AnalyticPricer::greeks: only vanilla calls and puts supported");
}

}  // namespace pricer
