#include "pricer/market_data.hpp"

#include <cmath>
#include <stdexcept>

namespace pricer {

MarketData::MarketData(double spot, double rate, double dividend, double volatility)
    : spot_(spot), rate_(rate), dividend_(dividend), volatility_(volatility) {
    if (!std::isfinite(spot) || !std::isfinite(rate) || !std::isfinite(dividend) ||
        !std::isfinite(volatility)) {
        throw std::invalid_argument("MarketData: all inputs must be finite");
    }
    if (spot <= 0.0) {
        throw std::invalid_argument("MarketData: spot must be > 0");
    }
    if (volatility < 0.0) {
        throw std::invalid_argument("MarketData: volatility must be >= 0");
    }
}

MarketData MarketData::withSpot(double spot) const {
    return MarketData(spot, rate_, dividend_, volatility_);
}

MarketData MarketData::withVolatility(double volatility) const {
    return MarketData(spot_, rate_, dividend_, volatility);
}

MarketData MarketData::withRate(double rate) const {
    return MarketData(spot_, rate, dividend_, volatility_);
}

}  // namespace pricer
