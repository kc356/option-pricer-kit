#pragma once

namespace OptionPricer
{

    /**
     * @brief Finite Difference method settings
     *
     * Controls grid parameters and solver settings for PDE-based option pricing
     */
    class FDSettings
    {
    public:
        /**
         * @brief Construct FD settings
         * @param S_max Upper bound of spot price domain (typically 3-5x strike)
         * @param N Number of space grid points (S direction)
         * @param M Number of time steps
         * @param theta Scheme parameter (0=explicit, 0.5=Crank-Nicolson, 1=implicit)
         * @param omega PSOR relaxation parameter (1.0-1.8 typical for American)
         * @param tol PSOR convergence tolerance
         * @param max_iter Maximum PSOR iterations per time step
         */
        FDSettings(double S_max = 300.0,
                   int N = 200,
                   int M = 200,
                   double theta = 0.5,
                   double omega = 1.2,
                   double tol = 1e-6,
                   int max_iter = 1000);

        // Getters
        double get_S_max() const { return S_max_; }
        int get_N() const { return N_; }
        int get_M() const { return M_; }
        double get_theta() const { return theta_; }
        double get_omega() const { return omega_; }
        double get_tol() const { return tol_; }
        int get_max_iter() const { return max_iter_; }

        // Setters
        void set_S_max(double S_max) { S_max_ = S_max; }
        void set_N(int N) { N_ = N; }
        void set_M(int M) { M_ = M; }
        void set_theta(double theta) { theta_ = theta; }
        void set_omega(double omega) { omega_ = omega; }
        void set_tol(double tol) { tol_ = tol; }
        void set_max_iter(int max_iter) { max_iter_ = max_iter; }

        // Computed properties
        double get_dS() const { return S_max_ / N_; }

    private:
        double S_max_; // Upper bound of S domain
        int N_;        // Number of space grid points
        int M_;        // Number of time steps
        double theta_; // Crank-Nicolson parameter (0.5 = CN)
        double omega_; // PSOR relaxation parameter
        double tol_;   // PSOR convergence tolerance
        int max_iter_; // Maximum PSOR iterations
    };

} // namespace OptionPricer
