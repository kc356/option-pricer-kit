#include "products/AsianOption.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace OptionPricer
{

    AsianOption::AsianOption(double K, double T, bool is_call, int n_steps)
        : K_(K), T_(T), is_call_(is_call), n_steps_(n_steps)
    {
        if (K <= 0.0)
        {
            throw std::invalid_argument("Strike must be positive");
        }
        if (T <= 0.0)
        {
            throw std::invalid_argument("Time to maturity must be positive");
        }
        if (n_steps <= 0)
        {
            throw std::invalid_argument("Number of steps must be positive");
        }
    }

    double AsianOption::payoff(const std::vector<double> &path) const
    {
        if (path.empty())
        {
            throw std::invalid_argument("Path cannot be empty");
        }

        // Arithmetic average
        double sum = std::accumulate(path.begin(), path.end(), 0.0);
        double avg = sum / path.size();

        if (is_call_)
        {
            return std::max(avg - K_, 0.0);
        }
        else
        {
            return std::max(K_ - avg, 0.0);
        }
    }

    double AsianOption::payoff_geometric(const std::vector<double> &path) const
    {
        if (path.empty())
        {
            throw std::invalid_argument("Path cannot be empty");
        }

        // Geometric average: exp(mean(log(S)))
        double log_sum = 0.0;
        for (double S : path)
        {
            if (S <= 0.0)
            {
                throw std::invalid_argument("Stock prices must be positive");
            }
            log_sum += std::log(S);
        }
        double geo_avg = std::exp(log_sum / path.size());

        if (is_call_)
        {
            return std::max(geo_avg - K_, 0.0);
        }
        else
        {
            return std::max(K_ - geo_avg, 0.0);
        }
    }

} // namespace OptionPricer
