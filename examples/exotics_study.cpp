// Bonus study: implied volatility solver, Asian options with a geometric control variate,
// and the discrete-monitoring bias of a barrier option (S=100, r=5%, q=2%, sigma=25%, T=1).
//
// Usage: exotics_study [output.csv]   (default: barrier_monitoring.csv)
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "pricer/black_scholes.hpp"
#include "pricer/exotic_formulas.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/path_dependent_pricer.hpp"

int main(int argc, char* argv[]) {
    using namespace pricer;

    const std::string path = argc > 1 ? argv[1] : "barrier_monitoring.csv";
    const MarketData market(100.0, 0.05, 0.02, 0.25);
    const black_scholes::Inputs atm{100.0, 100.0, 0.05, 0.02, 0.25, 1.0};
    constexpr double kMaturity = 1.0;
    constexpr std::size_t kPaths = 100000;
    constexpr std::uint64_t kSeed = 11;

    // --- 1. Implied volatility: price with a known vol, recover it ---------------------------
    std::cout << std::fixed << "1. Implied volatility (call, true sigma = 25%)\n"
              << std::setw(8) << "K" << std::setw(12) << "price" << std::setw(14) << "implied"
              << std::setw(8) << "iters" << '\n';
    for (double strike : {60.0, 80.0, 100.0, 120.0, 160.0}) {
        black_scholes::Inputs in = atm;
        in.strike = strike;
        const double price = black_scholes::vanillaPrice(OptionType::Call, in);
        const ImpliedVolResult r = impliedVolatility(OptionType::Call, price, in);
        std::cout << std::setprecision(0) << std::setw(8) << strike << std::setprecision(6)
                  << std::setw(12) << price << std::setprecision(10) << std::setw(14)
                  << r.volatility << std::setw(8) << r.iterations << '\n';
    }

    // --- 2. Arithmetic Asian call with and without the geometric control variate -------------
    constexpr std::size_t kFixings = 52;  // weekly fixings
    const PathDependentPricer weekly(kPaths, kFixings, kSeed);
    const ArithmeticAsianPayoff arithmetic(OptionType::Call, 100.0);
    const GeometricAsianPayoff geometric(OptionType::Call, 100.0);
    const double geometricExact = closed_form::geometricAsianPrice(OptionType::Call, atm, kFixings);
    const PricingResult plain = weekly.price(arithmetic, kMaturity, market);
    const PricingResult controlled =
        weekly.priceWithControlVariate(arithmetic, geometric, geometricExact, kMaturity, market);
    std::cout << std::setprecision(5)
              << "\n2. Arithmetic Asian call, K=100, 52 weekly fixings, N=" << kPaths
              << "\n   vanilla BS            " << black_scholes::vanillaPrice(OptionType::Call, atm)
              << "\n   geometric (exact)     " << geometricExact << "\n   arithmetic, plain MC  "
              << plain.price << "  (SE " << plain.standardError << ")"
              << "\n   arithmetic, geo CV    " << controlled.price << "  (SE "
              << controlled.standardError << ")\n   variance reduction    x" << std::setprecision(0)
              << (plain.standardError * plain.standardError) /
                     (controlled.standardError * controlled.standardError)
              << '\n';

    // --- 3. Down-and-out call: discrete monitoring bias and BGK correction -------------------
    const double barrier = 90.0;
    const BarrierPayoff downAndOut(OptionType::Call, 100.0, BarrierType::DownAndOut, barrier);
    const double continuous = closed_form::downAndOutCallPrice(atm, barrier);
    std::ofstream csv(path);
    if (!csv) {
        std::cerr << "Cannot open " << path << '\n';
        return 1;
    }
    csv << "monitoring_dates,mc_price,mc_std_error,continuous_price,bgk_price\n"
        << std::setprecision(10);
    std::cout << std::setprecision(5) << "\n3. Down-and-out call, K=100, H=90 (continuous price "
              << continuous << ")\n"
              << std::setw(8) << "M" << std::setw(12) << "MC" << std::setw(10) << "SE"
              << std::setw(12) << "BGK" << '\n';
    for (std::size_t m : {4U, 12U, 52U, 252U}) {
        const PricingResult mc =
            PathDependentPricer(kPaths, m, kSeed).price(downAndOut, kMaturity, market);
        const double bgk = closed_form::downAndOutCallPrice(
            atm, closed_form::continuityCorrectedBarrier(
                     barrier, 0.25, kMaturity / static_cast<double>(m), false));
        csv << m << ',' << mc.price << ',' << mc.standardError << ',' << continuous << ',' << bgk
            << '\n';
        std::cout << std::setw(8) << m << std::setw(12) << mc.price << std::setw(10)
                  << mc.standardError << std::setw(12) << bgk << '\n';
    }
    std::cout << "\nCSV written to " << path << '\n';
    return 0;
}
