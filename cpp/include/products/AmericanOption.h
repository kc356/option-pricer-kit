#pragma once

namespace OptionPricer
{

    /**
     * @brief American vanilla option (call or put) with early exercise
     */
    class AmericanOption
    {
    public:
        /**
         * @brief Construct American option
         * @param K Strike price
         * @param T Time to maturity (years)
         * @param is_call True for call, false for put
         */
        AmericanOption(double K, double T, bool is_call);

        // Getters
        double get_strike() const { return K_; }
        double get_maturity() const { return T_; }
        bool is_call() const { return is_call_; }

        /**
         * @brief Calculate payoff (same as European at maturity)
         * @param S Underlying price
         * @return Payoff value
         */
        double payoff(double S) const;

        /**
         * @brief Calculate immediate exercise value at any time
         * @param S Underlying price
         * @return Exercise value
         */
        double exercise_value(double S) const;

    private:
        double K_;     // Strike
        double T_;     // Time to maturity
        bool is_call_; // Call or put
    };

} // namespace OptionPricer
