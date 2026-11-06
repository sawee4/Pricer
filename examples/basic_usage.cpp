// Minimal tour of the building blocks: a market, a contract and its payoff.
// Pricing engines are added in the next steps of the roadmap.
#include <iostream>
#include <memory>

#include "pricer/market_data.hpp"
#include "pricer/option.hpp"

int main() {
    using namespace pricer;

    const MarketData market(/*spot=*/42.0, /*rate=*/0.10, /*dividend=*/0.0, /*vol=*/0.20);
    const Option call(std::make_unique<CallPayoff>(40.0), /*maturity=*/0.5);

    std::cout << "Spot           : " << market.spot() << '\n'
              << "Maturity (yrs) : " << call.maturity() << '\n'
              << "Intrinsic value: " << call.payoffAt(market.spot()) << '\n';
    return 0;
}
