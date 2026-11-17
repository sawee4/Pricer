# pricer — Option Pricing Library in Modern C++

[![CI](https://github.com/sawee4/pricer/actions/workflows/ci.yml/badge.svg)](https://github.com/sawee4/pricer/actions/workflows/ci.yml)

A small, tested option pricing library written in C++17, with no dependencies
other than GoogleTest (fetched automatically by CMake).

> Work in progress — see the roadmap below.

## Design

| Component     | Responsibility                                                       |
|---------------|----------------------------------------------------------------------|
| `Payoff`      | Terminal payoff `h(S_T)`: call, put, digital call/put (polymorphic, `clone()`) |
| `Option`      | Contract = payoff + maturity + exercise style (European / American)  |
| `MarketData`  | Black-Scholes inputs `S, r, q, σ`, validated at construction         |
| `Pricer`      | Pricing engine interface: analytic, Monte Carlo, binomial            |

## Build

Requirements: CMake ≥ 3.16 and a C++17 compiler (GCC, Clang, MSVC or MinGW).

```bash
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Roadmap

- [x] Project skeleton, CMake, CI (Linux GCC/Clang, Windows MSVC)
- [x] Black-Scholes closed form with continuous dividend yield (calls, puts, digitals)
- [x] Analytic Greeks (delta, gamma, vega, theta, rho), validated against finite differences
- [x] Monte Carlo pricer with standard error, 95% confidence interval and convergence study
- [x] Variance reduction: antithetic variates, control variate (with measured reduction factors)
- [ ] CRR binomial tree, European and American exercise
- [ ] Monte Carlo Greeks: bump-and-revalue with common random numbers, pathwise delta
- [ ] Implied volatility, exotic options (Asian, barrier)
