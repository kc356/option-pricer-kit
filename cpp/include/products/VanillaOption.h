#pragma once

namespace OptionPricer
{

    /**
     * @brief European vanilla option (call or put)
     */
    class VanillaOption
    {
    public:
        /**
         * @brief Construct vanilla option
         * @param K Strike price
         * @param T Time to maturity (years)
         * @param is_call True for call, false for put
         */
        VanillaOption(double K, double T, bool is_call);

        // Getters
        double get_strike() const { return K_; }
        double get_maturity() const { return T_; }
        bool is_call() const { return is_call_; }

        /**
         * @brief Calculate payoff at maturity
         * @param S_T Underlying price at maturity
         * @return Payoff value
         */
        double payoff(double S_T) const;

    private:
        double K_;     // Strike
        double T_;     // Time to maturity
        bool is_call_; // Call or put
    };

} // namespace OptionPricer
