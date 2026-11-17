// CRR binomial tree convergence study (S=100, K=105, r=5%, q=2%, sigma=25%, T=1).
//
// For n = 10..500 steps, records the European call price, its error against Black-Scholes,
// the error of the average of the n and n+1 trees (which damps the odd/even oscillation),
// and the American put price (no closed form; converges to a limit).
//
// Usage: binomial_convergence [output.csv]   (default: binomial_convergence.csv)
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "pricer/analytic_pricer.hpp"
#include "pricer/binomial_pricer.hpp"

int main(int argc, char* argv[]) {
    using namespace pricer;

    const std::string path = argc > 1 ? argv[1] : "binomial_convergence.csv";
    const MarketData market(100.0, 0.05, 0.02, 0.25);
    const Option euCall(std::make_unique<CallPayoff>(105.0), 1.0, ExerciseStyle::European);
    const Option usPut(std::make_unique<PutPayoff>(105.0), 1.0, ExerciseStyle::American);
    const Option euPut(std::make_unique<PutPayoff>(105.0), 1.0, ExerciseStyle::European);
    const double bs = AnalyticPricer().price(euCall, market).price;

    std::ofstream csv(path);
    if (!csv) {
        std::cerr << "Cannot open " << path << '\n';
        return 1;
    }
    csv << "n_steps,eu_call,bs_call,error,error_avg_n_n1,us_put,eu_put\n" << std::setprecision(10);

    std::cout << "European call BS price: " << std::fixed << std::setprecision(6) << bs << "\n\n"
              << std::setw(6) << "n" << std::setw(12) << "EU call" << std::setw(12) << "error"
              << std::setw(16) << "err avg(n,n+1)" << std::setw(12) << "US put" << std::setw(12)
              << "EU put" << '\n';

    for (std::size_t n = 10; n <= 500; ++n) {
        const double call = BinomialPricer(n).price(euCall, market).price;
        const double callNext = BinomialPricer(n + 1).price(euCall, market).price;
        const double us = BinomialPricer(n).price(usPut, market).price;
        const double eu = BinomialPricer(n).price(euPut, market).price;
        const double error = call - bs;
        const double errorAvg = 0.5 * (call + callNext) - bs;
        csv << n << ',' << call << ',' << bs << ',' << error << ',' << errorAvg << ',' << us << ','
            << eu << '\n';
        if (n == 10 || n == 11 || n == 50 || n == 51 || n == 100 || n == 101 || n == 500) {
            std::cout << std::setw(6) << n << std::setw(12) << call << std::setw(12) << error
                      << std::setw(16) << errorAvg << std::setw(12) << us << std::setw(12) << eu
                      << '\n';
        }
    }
    std::cout << "\nCSV written to " << path << '\n';
    return 0;
}
