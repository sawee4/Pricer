// Variance reduction study: plain Monte Carlo vs antithetic variates vs control variate
// (X = e^{-rT} S_T) for calls across moneyness, at a fixed number of payoff evaluations N.
//
// Reduction factor = Var(plain) / Var(method) = (SE_plain / SE_method)^2. Since every method
// uses the same N payoff evaluations, it is also the speed-up in paths needed for a given
// accuracy (the extra exp() per path of the control variate is negligible).
//
// Usage: variance_reduction_study [output.csv]   (default: variance_reduction.csv)
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "pricer/analytic_pricer.hpp"
#include "pricer/monte_carlo_pricer.hpp"

int main(int argc, char* argv[]) {
    using namespace pricer;

    const std::string path = argc > 1 ? argv[1] : "variance_reduction.csv";
    constexpr std::size_t kPaths = 100000;
    constexpr std::uint64_t kSeed = 2026;
    const MarketData market(100.0, 0.05, 0.0, 0.20);
    constexpr double kMaturity = 1.0;

    std::ofstream csv(path);
    if (!csv) {
        std::cerr << "Cannot open " << path << '\n';
        return 1;
    }
    csv << "strike,bs_price,se_plain,se_antithetic,se_control,factor_antithetic,factor_control\n";
    csv << std::setprecision(10);

    std::cout << "Calls, S=100, r=5%, sigma=20%, T=1, N=" << kPaths << " payoff evaluations\n\n"
              << std::setw(7) << "K" << std::setw(10) << "BS" << std::setw(11) << "SE plain"
              << std::setw(11) << "SE anti" << std::setw(11) << "SE ctrl" << std::setw(10)
              << "x anti" << std::setw(10) << "x ctrl" << '\n';

    for (double strike = 60.0; strike <= 160.0; strike += 10.0) {
        const Option call(std::make_unique<CallPayoff>(strike), kMaturity);
        const double bs = AnalyticPricer().price(call, market).price;
        const double sePlain = MonteCarloPricer(kPaths, kSeed).price(call, market).standardError;
        const double seAnti = MonteCarloPricer(kPaths, kSeed, VarianceReduction::Antithetic)
                                  .price(call, market).standardError;
        const double seCtrl = MonteCarloPricer(kPaths, kSeed, VarianceReduction::ControlVariate)
                                  .price(call, market).standardError;
        const double fAnti = (sePlain * sePlain) / (seAnti * seAnti);
        const double fCtrl = (sePlain * sePlain) / (seCtrl * seCtrl);

        csv << strike << ',' << bs << ',' << sePlain << ',' << seAnti << ',' << seCtrl << ','
            << fAnti << ',' << fCtrl << '\n';
        std::cout << std::fixed << std::setprecision(0) << std::setw(7) << strike
                  << std::setprecision(4) << std::setw(10) << bs << std::setw(11) << sePlain
                  << std::setw(11) << seAnti << std::setw(11) << seCtrl << std::setprecision(2)
                  << std::setw(10) << fAnti << std::setw(10) << fCtrl << '\n';
    }
    std::cout << "\nCSV written to " << path << '\n';
    return 0;
}
