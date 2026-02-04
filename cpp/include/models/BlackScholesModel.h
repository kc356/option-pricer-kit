#pragma once

namespace OptionPricer
{

    /**
     * @brief Black-Scholes model with constant volatility
     *
     * Represents the classic geometric Brownian motion model:
     * dS/S = (r - q)dt + sigma * dW
     */
    class BlackScholesModel
    {
    public:
        /**
         * @brief Construct Black-Scholes model
         * @param r Risk-free rate (continuously compounded)
         * @param q Dividend yield (continuously compounded)
         * @param sigma Volatility (annualized)
         */
        BlackScholesModel(double r, double q, double sigma);

        // Getters
        double get_r() const { return r_; }
        double get_q() const { return q_; }
        double get_sigma() const { return sigma_; }

        // Setters
        void set_r(double r) { r_ = r; }
        void set_q(double q) { q_ = q; }
        void set_sigma(double sigma) { sigma_ = sigma; }

    private:
        double r_;     // Risk-free rate
        double q_;     // Dividend yield
        double sigma_; // Volatility
    };

} // namespace OptionPricer
