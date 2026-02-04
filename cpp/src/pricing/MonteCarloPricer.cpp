#include "pricing/MonteCarloPricer.h"
#include "pricing/AnalyticPricer.h"
#include "utils/MathUtils.h"
#include <cmath>
#include <random>
#include <vector>
#include <numeric>

namespace OptionPricer
{

    // Helper function to generate a BS path using exact discretization
    static std::vector<double> generate_bs_path(
        double S0, double r, double q, double sigma, double T, int n_steps,
        std::mt19937 &rng, std::normal_distribution<> &dist)
    {

        std::vector<double> path;
        path.reserve(n_steps);

        double dt = T / n_steps;
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double vol = sigma * std::sqrt(dt);

        double S = S0;
        for (int i = 0; i < n_steps; ++i)
        {
            double Z = dist(rng);
            S = S * std::exp(drift + vol * Z);
            path.push_back(S);
        }

        return path;
    }

    MCResult mc_price(const VanillaOption &option,
                      const BlackScholesModel &model,
                      double S0,
                      const MCSettings &settings)
    {

        double K = option.get_strike();
        double T = option.get_maturity();
        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        int n_paths = settings.get_n_paths();
        int n_steps = settings.get_n_steps();
        bool antithetic = settings.use_antithetic();

        std::mt19937 rng(settings.get_seed());
        std::normal_distribution<> dist(0.0, 1.0);

        double dt = T / n_steps;
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double vol = sigma * std::sqrt(dt);
        double discount = std::exp(-r * T);

        std::vector<double> payoffs;
        payoffs.reserve(antithetic ? n_paths * 2 : n_paths);

        for (int i = 0; i < n_paths; ++i)
        {
            // Generate path
            double S = S0;
            for (int j = 0; j < n_steps; ++j)
            {
                double Z = dist(rng);
                S = S * std::exp(drift + vol * Z);
            }
            payoffs.push_back(option.payoff(S));

            // Antithetic path
            if (antithetic)
            {
                S = S0;
                for (int j = 0; j < n_steps; ++j)
                {
                    double Z = dist(rng);
                    S = S * std::exp(drift - vol * Z); // Note the minus sign
                }
                payoffs.push_back(option.payoff(S));
            }
        }

        // Calculate statistics
        double sum = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
        double mean = sum / payoffs.size();

        double sq_sum = 0.0;
        for (double p : payoffs)
        {
            sq_sum += (p - mean) * (p - mean);
        }
        double variance = sq_sum / (payoffs.size() - 1);
        double std_error = std::sqrt(variance / payoffs.size());

        MCResult result;
        result.price = discount * mean;
        result.std_error = discount * std_error;
        result.ci_lower = result.price - 1.96 * result.std_error;
        result.ci_upper = result.price + 1.96 * result.std_error;

        return result;
    }

    MCResult mc_price(const AsianOption &option,
                      const BlackScholesModel &model,
                      double S0,
                      const MCSettings &settings)
    {

        double K = option.get_strike();
        double T = option.get_maturity();
        int n_monitoring = option.get_n_steps();
        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        int n_paths = settings.get_n_paths();
        int n_steps = settings.get_n_steps();
        bool antithetic = settings.use_antithetic();

        std::mt19937 rng(settings.get_seed());
        std::normal_distribution<> dist(0.0, 1.0);

        double discount = std::exp(-r * T);

        std::vector<double> payoffs;
        payoffs.reserve(antithetic ? n_paths * 2 : n_paths);

        for (int i = 0; i < n_paths; ++i)
        {
            // Generate path and extract monitoring points
            auto full_path = generate_bs_path(S0, r, q, sigma, T, n_monitoring, rng, dist);
            payoffs.push_back(option.payoff(full_path));

            // Antithetic path
            if (antithetic)
            {
                // Re-seed to get same randoms but negated
                std::mt19937 rng_anti(settings.get_seed() + i + 1);
                std::normal_distribution<> dist_anti(0.0, 1.0);

                double dt = T / n_monitoring;
                double drift = (r - q - 0.5 * sigma * sigma) * dt;
                double vol = sigma * std::sqrt(dt);

                std::vector<double> anti_path;
                anti_path.reserve(n_monitoring);
                double S = S0;
                for (int j = 0; j < n_monitoring; ++j)
                {
                    double Z = dist(rng); // Use original Z
                    S = S0;               // Reset
                    // Recalculate with accumulated -Z
                    for (int k = 0; k <= j; ++k)
                    {
                        double Z_k = dist(rng);
                        if (k == j)
                        {
                            S = S * std::exp(drift - vol * Z_k);
                        }
                    }
                }

                // Simpler approach: just negate the increments
                S = S0;
                std::vector<double> zs;
                for (int j = 0; j < n_monitoring; ++j)
                {
                    double Z = dist(rng);
                    zs.push_back(Z);
                }

                std::vector<double> anti_path2;
                S = S0;
                for (int j = 0; j < n_monitoring; ++j)
                {
                    S = S * std::exp(drift - vol * zs[j]);
                    anti_path2.push_back(S);
                }
                payoffs.push_back(option.payoff(anti_path2));
            }
        }

        // Calculate statistics
        double sum = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
        double mean = sum / payoffs.size();

        double sq_sum = 0.0;
        for (double p : payoffs)
        {
            sq_sum += (p - mean) * (p - mean);
        }
        double variance = sq_sum / (payoffs.size() - 1);
        double std_error = std::sqrt(variance / payoffs.size());

        MCResult result;
        result.price = discount * mean;
        result.std_error = discount * std_error;
        result.ci_lower = result.price - 1.96 * result.std_error;
        result.ci_upper = result.price + 1.96 * result.std_error;

        return result;
    }

    MCResult mc_price_asian_cv(const AsianOption &option,
                               const BlackScholesModel &model,
                               double S0,
                               const MCSettings &settings,
                               bool use_cv)
    {

        if (!use_cv)
        {
            return mc_price(option, model, S0, settings);
        }

        double K = option.get_strike();
        double T = option.get_maturity();
        int n_monitoring = option.get_n_steps();
        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        int n_paths = settings.get_n_paths();
        bool antithetic = settings.use_antithetic();

        std::mt19937 rng(settings.get_seed());
        std::normal_distribution<> dist(0.0, 1.0);

        double discount = std::exp(-r * T);

        // Calculate control variate expected value (geometric Asian analytic price)
        double cv_expected = price_geometric_asian_analytic(K, T, option.is_call(),
                                                            n_monitoring, model, S0) /
                             discount;

        std::vector<double> arithmetic_payoffs;
        std::vector<double> geometric_payoffs;

        int total_paths = antithetic ? n_paths * 2 : n_paths;
        arithmetic_payoffs.reserve(total_paths);
        geometric_payoffs.reserve(total_paths);

        for (int i = 0; i < n_paths; ++i)
        {
            // Generate path
            auto path = generate_bs_path(S0, r, q, sigma, T, n_monitoring, rng, dist);
            arithmetic_payoffs.push_back(option.payoff(path));
            geometric_payoffs.push_back(option.payoff_geometric(path));

            // Antithetic path
            if (antithetic)
            {
                double dt = T / n_monitoring;
                double drift = (r - q - 0.5 * sigma * sigma) * dt;
                double vol = sigma * std::sqrt(dt);

                std::vector<double> anti_path;
                anti_path.reserve(n_monitoring);
                double S = S0;
                for (size_t j = 0; j < path.size(); ++j)
                {
                    double Z = dist(rng);
                    S = S * std::exp(drift - vol * Z);
                    anti_path.push_back(S);
                }
                arithmetic_payoffs.push_back(option.payoff(anti_path));
                geometric_payoffs.push_back(option.payoff_geometric(anti_path));
            }
        }

        // Calculate covariance and optimal beta
        double mean_arith = std::accumulate(arithmetic_payoffs.begin(), arithmetic_payoffs.end(), 0.0) / arithmetic_payoffs.size();
        double mean_geom = std::accumulate(geometric_payoffs.begin(), geometric_payoffs.end(), 0.0) / geometric_payoffs.size();

        double cov = 0.0;
        double var_geom = 0.0;
        for (size_t i = 0; i < arithmetic_payoffs.size(); ++i)
        {
            cov += (arithmetic_payoffs[i] - mean_arith) * (geometric_payoffs[i] - mean_geom);
            var_geom += (geometric_payoffs[i] - mean_geom) * (geometric_payoffs[i] - mean_geom);
        }
        cov /= (arithmetic_payoffs.size() - 1);
        var_geom /= (arithmetic_payoffs.size() - 1);

        double beta = (var_geom > 0) ? cov / var_geom : 0.0;

        // Apply control variate
        std::vector<double> cv_payoffs;
        cv_payoffs.reserve(arithmetic_payoffs.size());
        for (size_t i = 0; i < arithmetic_payoffs.size(); ++i)
        {
            double cv_payoff = arithmetic_payoffs[i] - beta * (geometric_payoffs[i] - cv_expected);
            cv_payoffs.push_back(cv_payoff);
        }

        // Calculate statistics
        double sum = std::accumulate(cv_payoffs.begin(), cv_payoffs.end(), 0.0);
        double mean = sum / cv_payoffs.size();

        double sq_sum = 0.0;
        for (double p : cv_payoffs)
        {
            sq_sum += (p - mean) * (p - mean);
        }
        double variance = sq_sum / (cv_payoffs.size() - 1);
        double std_error = std::sqrt(variance / cv_payoffs.size());

        MCResult result;
        result.price = discount * mean;
        result.std_error = discount * std_error;
        result.ci_lower = result.price - 1.96 * result.std_error;
        result.ci_upper = result.price + 1.96 * result.std_error;

        return result;
    }

} // namespace OptionPricer
