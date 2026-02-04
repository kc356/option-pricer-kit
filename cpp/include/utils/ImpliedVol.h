#pragma once

#include <stdexcept>

namespace OptionPricer
{

    /**
     * @brief Calculate implied volatility using bisection method
     *
     * Inverts the Black-Scholes formula to find the volatility that
     * produces the given market price
     *
     * @param market_price Observed option price in market
     * @param S0 Current underlying price
     * @param K Strike price
     * @param T Time to maturity (years)
     * @param r Risk-free rate
     * @param q Dividend yield
     * @param is_call True for call, false for put
     * @param tol Convergence tolerance (default 1e-6)
     * @param max_iter Maximum iterations (default 100)
     * @return Implied volatility
     * @throws std::runtime_error if convergence fails or price is invalid
     */
    double implied_vol_bisect(double market_price,
                              double S0, double K, double T,
                              double r, double q, bool is_call,
                              double tol = 1e-6,
                              int max_iter = 100);

} // namespace OptionPricer
