#pragma once

#include <vector>

namespace OptionPricer
{

    /**
     * @brief Asian option with arithmetic average (discrete monitoring)
     */
    class AsianOption
    {
    public:
        /**
         * @brief Construct Asian option
         * @param K Strike price
         * @param T Time to maturity (years)
         * @param is_call True for call, false for put
         * @param n_steps Number of monitoring points (equally spaced)
         */
        AsianOption(double K, double T, bool is_call, int n_steps);

        // Getters
        double get_strike() const { return K_; }
        double get_maturity() const { return T_; }
        bool is_call() const { return is_call_; }
        int get_n_steps() const { return n_steps_; }

        /**
         * @brief Calculate payoff based on arithmetic average
         * @param path Vector of underlying prices at monitoring points
         * @return Payoff value
         */
        double payoff(const std::vector<double> &path) const;

        /**
         * @brief Calculate payoff based on geometric average
         * @param path Vector of underlying prices at monitoring points
         * @return Payoff value (used for control variate)
         */
        double payoff_geometric(const std::vector<double> &path) const;

    private:
        double K_;     // Strike
        double T_;     // Time to maturity
        bool is_call_; // Call or put
        int n_steps_;  // Number of monitoring points
    };

} // namespace OptionPricer
