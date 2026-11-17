#include "pricer/exotic_formulas.hpp"

#include <cmath>
#include <stdexcept>

#include "pricer/normal.hpp"

namespace pricer::closed_form {

double geometricAsianPrice(OptionType type, const black_scholes::Inputs& in,
                           std::size_t numFixings) {
    if (numFixings == 0 || !(in.maturity > 0.0) || !(in.volatility > 0.0)) {
        throw std::invalid_argument("geometricAsianPrice: need M >= 1, T > 0, sigma > 0");
    }
    const double m = static_cast<double>(numFixings);
    const double sigma2 = in.volatility * in.volatility;
    const double mu = std::log(in.spot) +
                      (in.rate - in.dividend - 0.5 * sigma2) * in.maturity * (m + 1.0) / (2.0 * m);
    const double variance = sigma2 * in.maturity * (m + 1.0) * (2.0 * m + 1.0) / (6.0 * m * m);
    const double s = std::sqrt(variance);
    const double d1 = (mu - std::log(in.strike) + variance) / s;
    const double d2 = d1 - s;
    const double expectedG = std::exp(mu + 0.5 * variance);  // E[G]
    const double df = std::exp(-in.rate * in.maturity);
    if (type == OptionType::Call) {
        return df * (expectedG * math::normCdf(d1) - in.strike * math::normCdf(d2));
    }
    return df * (in.strike * math::normCdf(-d2) - expectedG * math::normCdf(-d1));
}

double downAndInCallPrice(const black_scholes::Inputs& in, double barrier) {
    if (!(barrier > 0.0) || barrier > in.strike || !(barrier < in.spot) || !(in.maturity > 0.0) ||
        !(in.volatility > 0.0)) {
        throw std::invalid_argument("downAndInCallPrice: need 0 < H <= K, H < S, T > 0, sigma > 0");
    }
    const double sigma2 = in.volatility * in.volatility;
    const double volSqrtT = in.volatility * std::sqrt(in.maturity);
    const double lambda = (in.rate - in.dividend + 0.5 * sigma2) / sigma2;
    const double ratio = barrier / in.spot;
    const double y = std::log(barrier * barrier / (in.spot * in.strike)) / volSqrtT +
                     lambda * volSqrtT;
    return in.spot * std::exp(-in.dividend * in.maturity) * std::pow(ratio, 2.0 * lambda) *
               math::normCdf(y) -
           in.strike * std::exp(-in.rate * in.maturity) * std::pow(ratio, 2.0 * lambda - 2.0) *
               math::normCdf(y - volSqrtT);
}

double downAndOutCallPrice(const black_scholes::Inputs& in, double barrier) {
    return black_scholes::vanillaPrice(OptionType::Call, in) - downAndInCallPrice(in, barrier);
}

double continuityCorrectedBarrier(double barrier, double volatility, double dt, bool upBarrier) {
    const double shift = kBgkBeta * volatility * std::sqrt(dt);
    return barrier * std::exp(upBarrier ? shift : -shift);
}

}  // namespace pricer::closed_form
