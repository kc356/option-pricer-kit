#include "pricing/AnalyticPricer.h"
#include "utils/MathUtils.h"
#include <cmath>
#include <stdexcept>

namespace OptionPricer
{

    double price_bs_analytic(const VanillaOption &option,
                             const BlackScholesModel &model,
                             double S0)
    {
        if (S0 <= 0.0)
        {
            throw std::invalid_argument("Initial price must be positive");
        }

        double K = option.get_strike();
        double T = option.get_maturity();
        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        // Black-Scholes d1 and d2
        double sigma_sqrt_T = sigma * std::sqrt(T);
        double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / sigma_sqrt_T;
        double d2 = d1 - sigma_sqrt_T;

        double price;
        if (option.is_call())
        {
            price = S0 * std::exp(-q * T) * MathUtils::norm_cdf(d1) - K * std::exp(-r * T) * MathUtils::norm_cdf(d2);
        }
        else
        {
            price = K * std::exp(-r * T) * MathUtils::norm_cdf(-d2) - S0 * std::exp(-q * T) * MathUtils::norm_cdf(-d1);
        }

        return price;
    }

    double price_geometric_asian_analytic(double K, double T, bool is_call,
                                          int n_steps,
                                          const BlackScholesModel &model,
                                          double S0)
    {
        if (S0 <= 0.0)
        {
            throw std::invalid_argument("Initial price must be positive");
        }
        if (n_steps <= 0)
        {
            throw std::invalid_argument("Number of steps must be positive");
        }

        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        // Adjusted parameters for geometric Asian
        // The geometric average of lognormal variables is lognormal with adjusted parameters
        double sigma_adj = sigma * std::sqrt((n_steps + 1.0) * (2.0 * n_steps + 1.0) / (6.0 * n_steps * n_steps));
        double mu = (r - q - 0.5 * sigma * sigma) * (n_steps + 1.0) / (2.0 * n_steps);
        double r_adj = mu + 0.5 * sigma_adj * sigma_adj;

        // Use BS formula with adjusted parameters
        double sigma_sqrt_T = sigma_adj * std::sqrt(T);
        double d1 = (std::log(S0 / K) + (r_adj + 0.5 * sigma_adj * sigma_adj) * T) / sigma_sqrt_T;
        double d2 = d1 - sigma_sqrt_T;

        double price;
        if (is_call)
        {
            price = S0 * std::exp((mu - r) * T) * MathUtils::norm_cdf(d1) - K * std::exp(-r * T) * MathUtils::norm_cdf(d2);
        }
        else
        {
            price = K * std::exp(-r * T) * MathUtils::norm_cdf(-d2) - S0 * std::exp((mu - r) * T) * MathUtils::norm_cdf(-d1);
        }

        return price;
    }

} // namespace OptionPricer
