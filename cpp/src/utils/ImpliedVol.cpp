#include "utils/ImpliedVol.h"
#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"
#include "pricing/AnalyticPricer.h"
#include <cmath>
#include <limits>

namespace OptionPricer
{

    double implied_vol_bisect(double market_price,
                              double S0, double K, double T,
                              double r, double q, bool is_call,
                              double tol,
                              int max_iter)
    {

        // Validate inputs
        if (market_price <= 0.0)
        {
            throw std::runtime_error("Market price must be positive");
        }
        if (S0 <= 0.0 || K <= 0.0 || T <= 0.0)
        {
            throw std::runtime_error("S0, K, and T must be positive");
        }

        // Check intrinsic value
        double intrinsic = is_call ? std::max(S0 * std::exp(-q * T) - K * std::exp(-r * T), 0.0)
                                   : std::max(K * std::exp(-r * T) - S0 * std::exp(-q * T), 0.0);
        if (market_price < intrinsic - tol)
        {
            throw std::runtime_error("Market price below intrinsic value");
        }

        // Initial bracket
        double sigma_low = 1e-6;
        double sigma_high = 5.0;

        VanillaOption option(K, T, is_call);

        // Check if bracket is valid
        BlackScholesModel model_low(r, q, sigma_low);
        BlackScholesModel model_high(r, q, sigma_high);

        double price_low = price_bs_analytic(option, model_low, S0);
        double price_high = price_bs_analytic(option, model_high, S0);

        if (market_price < price_low || market_price > price_high)
        {
            // Market price outside reasonable volatility range
            if (market_price < price_low)
            {
                return sigma_low;
            }
            else
            {
                return sigma_high;
            }
        }

        // Bisection method
        for (int iter = 0; iter < max_iter; ++iter)
        {
            double sigma_mid = 0.5 * (sigma_low + sigma_high);
            BlackScholesModel model_mid(r, q, sigma_mid);
            double price_mid = price_bs_analytic(option, model_mid, S0);

            if (std::abs(price_mid - market_price) < tol)
            {
                return sigma_mid;
            }

            if (price_mid < market_price)
            {
                sigma_low = sigma_mid;
            }
            else
            {
                sigma_high = sigma_mid;
            }

            if (sigma_high - sigma_low < tol)
            {
                return 0.5 * (sigma_low + sigma_high);
            }
        }

        // Max iterations reached
        return 0.5 * (sigma_low + sigma_high);
    }

} // namespace OptionPricer
