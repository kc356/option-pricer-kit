#pragma once

#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"

namespace OptionPricer
{

    /**
     * @brief Black-Scholes analytic pricing for European vanilla options
     * @param option Vanilla option to price
     * @param model Black-Scholes model parameters
     * @param S0 Initial underlying price
     * @return Option price
     */
    double price_bs_analytic(const VanillaOption &option,
                             const BlackScholesModel &model,
                             double S0);

    /**
     * @brief Closed-form price for geometric Asian option under BS
     * @param K Strike price
     * @param T Time to maturity
     * @param is_call True for call, false for put
     * @param n_steps Number of monitoring points
     * @param model Black-Scholes model parameters
     * @param S0 Initial underlying price
     * @return Geometric Asian option price
     */
    double price_geometric_asian_analytic(double K, double T, bool is_call,
                                          int n_steps,
                                          const BlackScholesModel &model,
                                          double S0);

} // namespace OptionPricer
