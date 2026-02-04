#pragma once

namespace OptionPricer
{

    /**
     * @brief Barrier option types
     */
    enum class BarrierType
    {
        UpAndOut,   // Knock-out when S >= B
        DownAndOut, // Knock-out when S <= B
        UpAndIn,    // Knock-in when S >= B
        DownAndIn   // Knock-in when S <= B
    };

    /**
     * @brief European barrier option (call or put)
     *
     * Continuous monitoring barrier (for FD pricing).
     * - Knock-out: Option becomes worthless if barrier is breached
     * - Knock-in: Option only activates if barrier is breached
     */
    class BarrierOption
    {
    public:
        /**
         * @brief Construct barrier option
         * @param K Strike price
         * @param T Time to maturity (years)
         * @param is_call True for call, false for put
         * @param barrier_type Type of barrier (Up/Down, In/Out)
         * @param B Barrier level
         * @param rebate Rebate paid on knock-out (default 0)
         */
        BarrierOption(double K, double T, bool is_call,
                      BarrierType barrier_type, double B,
                      double rebate = 0.0);

        // Getters
        double get_strike() const { return K_; }
        double get_maturity() const { return T_; }
        bool is_call() const { return is_call_; }
        BarrierType get_barrier_type() const { return barrier_type_; }
        double get_barrier() const { return B_; }
        double get_rebate() const { return rebate_; }

        /**
         * @brief Check if barrier type is knock-out
         */
        bool is_knock_out() const
        {
            return barrier_type_ == BarrierType::UpAndOut ||
                   barrier_type_ == BarrierType::DownAndOut;
        }

        /**
         * @brief Check if barrier type is knock-in
         */
        bool is_knock_in() const
        {
            return barrier_type_ == BarrierType::UpAndIn ||
                   barrier_type_ == BarrierType::DownAndIn;
        }

        /**
         * @brief Check if barrier is up-type
         */
        bool is_up() const
        {
            return barrier_type_ == BarrierType::UpAndOut ||
                   barrier_type_ == BarrierType::UpAndIn;
        }

        /**
         * @brief Check if barrier is down-type
         */
        bool is_down() const
        {
            return barrier_type_ == BarrierType::DownAndOut ||
                   barrier_type_ == BarrierType::DownAndIn;
        }

        /**
         * @brief Calculate payoff at maturity (assuming not knocked out)
         * @param S_T Underlying price at maturity
         * @return Payoff value
         */
        double payoff(double S_T) const;

        /**
         * @brief Check if barrier has been breached
         * @param S Current underlying price
         * @return True if barrier breached
         */
        bool is_breached(double S) const;

    private:
        double K_;                 // Strike
        double T_;                 // Time to maturity
        bool is_call_;             // Call or put
        BarrierType barrier_type_; // Barrier type
        double B_;                 // Barrier level
        double rebate_;            // Rebate on knock-out
    };

    /**
     * @brief Convert barrier type to string
     */
    const char *barrier_type_to_string(BarrierType type);

} // namespace OptionPricer
