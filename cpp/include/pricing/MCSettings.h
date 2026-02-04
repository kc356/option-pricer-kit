#pragma once

namespace OptionPricer
{

    /**
     * @brief Monte Carlo simulation settings
     */
    class MCSettings
    {
    public:
        /**
         * @brief Construct MC settings
         * @param n_paths Number of simulation paths
         * @param n_steps Number of time steps per path
         * @param seed Random seed (0 for random)
         * @param antithetic Use antithetic variates for variance reduction
         */
        MCSettings(int n_paths = 100000, int n_steps = 252,
                   unsigned int seed = 42, bool antithetic = true);

        // Getters
        int get_n_paths() const { return n_paths_; }
        int get_n_steps() const { return n_steps_; }
        unsigned int get_seed() const { return seed_; }
        bool use_antithetic() const { return antithetic_; }

        // Setters
        void set_n_paths(int n) { n_paths_ = n; }
        void set_n_steps(int n) { n_steps_ = n; }
        void set_seed(unsigned int s) { seed_ = s; }
        void set_antithetic(bool a) { antithetic_ = a; }

    private:
        int n_paths_;       // Number of paths
        int n_steps_;       // Time steps per path
        unsigned int seed_; // Random seed
        bool antithetic_;   // Antithetic variates flag
    };

} // namespace OptionPricer
