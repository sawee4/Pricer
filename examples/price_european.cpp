// Prices the textbook example from Hull (Example 15.6) with the analytic pricer:
// S = 42, K = 40, r = 10%, sigma = 20%, T = 0.5  ->  call ~ 4.76, put ~ 0.81.
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

#include "pricer/analytic_pricer.hpp"

int main() {
    using namespace pricer;

    const MarketData market(/*spot=*/42.0, /*rate=*/0.10, /*dividend=*/0.0, /*vol=*/0.20);
    const Option call(std::make_unique<CallPayoff>(40.0), /*maturity=*/0.5);
    const Option put(std::make_unique<PutPayoff>(40.0), /*maturity=*/0.5);

    const AnalyticPricer pricer;
    const double c = pricer.price(call, market).price;
    const double p = pricer.price(put, market).price;

    std::cout << std::fixed << std::setprecision(6)
              << "European call : " << c << '\n'
              << "European put  : " << p << '\n'
              << "C - P         : " << c - p << '\n'
              << "S - K e^{-rT} : " << 42.0 - 40.0 * std::exp(-0.10 * 0.5) << '\n';

    const Greeks g = pricer.greeks(call, market);
    std::cout << "\nCall Greeks (vega/rho per 1.0, theta per year)\n"
              << "  delta : " << g.delta << '\n'
              << "  gamma : " << g.gamma << '\n'
              << "  vega  : " << g.vega << '\n'
              << "  theta : " << g.theta << '\n'
              << "  rho   : " << g.rho << '\n';
    return 0;
}
