#include "pricing/FDPricer.h"
#include "pricing/AnalyticPricer.h"
#include "utils/FDUtils.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace OptionPricer
{

    // ========================================================================
    // Internal helper functions
    // ========================================================================

    namespace
    {

        /**
         * @brief Build tridiagonal coefficients for BS PDE (Crank-Nicolson)
         *
         * The BS PDE: dV/dt + 0.5*sigma^2*S^2*d2V/dS2 + (r-q)*S*dV/dS - r*V = 0
         *
         * Discretized on uniform S grid with central differences.
         *
         * @param S Grid values
         * @param r Risk-free rate
         * @param q Dividend yield
         * @param sigma Volatility
         * @param dt Time step
         * @param theta Scheme parameter (0.5 = Crank-Nicolson)
         * @param a_impl Sub-diagonal (implicit part, LHS)
         * @param b_impl Main diagonal (implicit, LHS)
         * @param c_impl Super-diagonal (implicit, LHS)
         * @param a_expl Sub-diagonal (explicit part, RHS)
         * @param b_expl Main diagonal (explicit, RHS)
         * @param c_expl Super-diagonal (explicit, RHS)
         */
        void build_tridiagonal_coeffs(const std::vector<double> &S,
                                      double r, double q, double sigma,
                                      double dt, double theta,
                                      std::vector<double> &a_impl,
                                      std::vector<double> &b_impl,
                                      std::vector<double> &c_impl,
                                      std::vector<double> &a_expl,
                                      std::vector<double> &b_expl,
                                      std::vector<double> &c_expl)
        {
            int N = static_cast<int>(S.size());
            double dS = S[1] - S[0]; // Assume uniform grid

            a_impl.resize(N - 1);
            b_impl.resize(N);
            c_impl.resize(N - 1);
            a_expl.resize(N - 1);
            b_expl.resize(N);
            c_expl.resize(N - 1);

            for (int i = 1; i < N - 1; ++i)
            {
                double Si = S[i];
                double sigma2 = sigma * sigma;

                // Coefficients for -0.5*sigma^2*S^2*d2V/dS2 - (r-q)*S*dV/dS + r*V
                // Using central differences
                double alpha = 0.5 * dt * (sigma2 * Si * Si / (dS * dS) - (r - q) * Si / dS);
                double beta = dt * (sigma2 * Si * Si / (dS * dS) + r);
                double gamma = 0.5 * dt * (sigma2 * Si * Si / (dS * dS) + (r - q) * Si / dS);

                // Implicit part (time n+1): (I + theta * A)
                a_impl[i - 1] = -theta * alpha;
                b_impl[i] = 1.0 + theta * beta;
                c_impl[i] = -theta * gamma;

                // Explicit part (time n): (I - (1-theta) * A)
                a_expl[i - 1] = (1.0 - theta) * alpha;
                b_expl[i] = 1.0 - (1.0 - theta) * beta;
                c_expl[i] = (1.0 - theta) * gamma;
            }

            // Boundary conditions will be handled separately
            // Initialize boundary coefficients as identity
            b_impl[0] = 1.0;
            c_impl[0] = 0.0;
            b_impl[N - 1] = 1.0;
            a_impl[N - 2] = 0.0;

            b_expl[0] = 1.0;
            c_expl[0] = 0.0;
            b_expl[N - 1] = 1.0;
            a_expl[N - 2] = 0.0;
        }

        /**
         * @brief Apply matrix-vector product for explicit part
         */
        void apply_explicit(const std::vector<double> &a,
                            const std::vector<double> &b,
                            const std::vector<double> &c,
                            const std::vector<double> &V,
                            std::vector<double> &result)
        {
            int N = static_cast<int>(V.size());
            result.resize(N);

            result[0] = b[0] * V[0] + c[0] * V[1];
            for (int i = 1; i < N - 1; ++i)
            {
                result[i] = a[i - 1] * V[i - 1] + b[i] * V[i] + c[i] * V[i + 1];
            }
            result[N - 1] = a[N - 2] * V[N - 2] + b[N - 1] * V[N - 1];
        }

    } // anonymous namespace

    // ========================================================================
    // European Vanilla FD
    // ========================================================================

    FDResult fd_price_european(const VanillaOption &option,
                               const BlackScholesModel &model,
                               double S0,
                               const FDSettings &settings)
    {
        double K = option.get_strike();
        double T = option.get_maturity();
        bool is_call = option.is_call();

        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        double S_max = settings.get_S_max();
        int N = settings.get_N();
        int M = settings.get_M();
        double theta = settings.get_theta();

        // Build grid
        std::vector<double> S = uniform_grid(S_max, N);
        double dS = S[1] - S[0];
        double dt = T / M;

        // Initialize option values at maturity
        std::vector<double> V(N);
        std::vector<double> payoff_vec(N);
        for (int i = 0; i < N; ++i)
        {
            V[i] = option.payoff(S[i]);
            payoff_vec[i] = V[i];
        }

        // Build coefficient matrices
        std::vector<double> a_impl, b_impl, c_impl;
        std::vector<double> a_expl, b_expl, c_expl;
        build_tridiagonal_coeffs(S, r, q, sigma, dt, theta,
                                 a_impl, b_impl, c_impl,
                                 a_expl, b_expl, c_expl);

        // Time stepping (backward from T to 0)
        std::vector<double> rhs(N);
        std::vector<double> V_new(N);

        for (int m = M - 1; m >= 0; --m)
        {
            double t = m * dt; // Current time (going backward)

            // Apply explicit operator to get RHS
            apply_explicit(a_expl, b_expl, c_expl, V, rhs);

            // Apply boundary conditions
            if (is_call)
            {
                // Call: V(0,t) = 0
                rhs[0] = 0.0;
                // V(S_max,t) ~ S_max*exp(-q*t) - K*exp(-r*t)
                rhs[N - 1] = S_max * std::exp(-q * t) - K * std::exp(-r * t);
                if (rhs[N - 1] < 0.0)
                    rhs[N - 1] = 0.0;
            }
            else
            {
                // Put: V(0,t) ~ K*exp(-r*t)
                rhs[0] = K * std::exp(-r * t);
                // V(S_max,t) = 0
                rhs[N - 1] = 0.0;
            }

            // Solve implicit system
            solve_tridiagonal(a_impl, b_impl, c_impl, rhs, V_new);

            // Enforce boundary conditions on solution
            if (is_call)
            {
                V_new[0] = 0.0;
                V_new[N - 1] = std::max(S_max * std::exp(-q * t) - K * std::exp(-r * t), 0.0);
            }
            else
            {
                V_new[0] = K * std::exp(-r * t);
                V_new[N - 1] = 0.0;
            }

            V = V_new;
        }

        // Interpolate to get V(S0)
        double price = linear_interp(S, V, S0);

        FDResult result;
        result.price = price;
        result.grid_S = S;
        result.grid_V = V;
        result.grid_payoff = payoff_vec;

        return result;
    }

    // ========================================================================
    // American Vanilla FD with PSOR
    // ========================================================================

    FDResult fd_price_american(const AmericanOption &option,
                               const BlackScholesModel &model,
                               double S0,
                               const FDSettings &settings)
    {
        double K = option.get_strike();
        double T = option.get_maturity();
        bool is_call = option.is_call();

        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        double S_max = settings.get_S_max();
        int N = settings.get_N();
        int M = settings.get_M();
        double theta = settings.get_theta();
        double omega = settings.get_omega();
        double tol = settings.get_tol();
        int max_iter = settings.get_max_iter();

        // Build grid
        std::vector<double> S = uniform_grid(S_max, N);
        double dS = S[1] - S[0];
        double dt = T / M;

        // Initialize option values at maturity
        std::vector<double> V(N);
        std::vector<double> payoff_vec(N);
        for (int i = 0; i < N; ++i)
        {
            V[i] = option.payoff(S[i]);
            payoff_vec[i] = V[i];
        }

        // Build coefficient matrices
        std::vector<double> a_impl, b_impl, c_impl;
        std::vector<double> a_expl, b_expl, c_expl;
        build_tridiagonal_coeffs(S, r, q, sigma, dt, theta,
                                 a_impl, b_impl, c_impl,
                                 a_expl, b_expl, c_expl);

        // Exercise value (lower bound for American)
        std::vector<double> exercise(N);
        for (int i = 0; i < N; ++i)
        {
            exercise[i] = option.exercise_value(S[i]);
        }

        // Time stepping (backward from T to 0)
        std::vector<double> rhs(N);

        for (int m = M - 1; m >= 0; --m)
        {
            double t = m * dt;

            // Apply explicit operator to get RHS
            apply_explicit(a_expl, b_expl, c_expl, V, rhs);

            // Apply boundary conditions
            if (is_call)
            {
                rhs[0] = 0.0;
                exercise[0] = 0.0;
                rhs[N - 1] = std::max(S_max - K, 0.0);
                exercise[N - 1] = std::max(S_max - K, 0.0);
            }
            else
            {
                rhs[0] = K; // For American put, at S=0, value = K
                exercise[0] = K;
                rhs[N - 1] = 0.0;
                exercise[N - 1] = 0.0;
            }

            // Solve with PSOR (enforces V >= exercise)
            solve_psor(a_impl, b_impl, c_impl, rhs, exercise, V, omega, tol, max_iter);

            // Enforce boundary values
            if (is_call)
            {
                V[0] = 0.0;
                V[N - 1] = std::max(S_max - K, 0.0);
            }
            else
            {
                V[0] = K;
                V[N - 1] = 0.0;
            }
        }

        // Interpolate to get V(S0)
        double price = linear_interp(S, V, S0);

        FDResult result;
        result.price = price;
        result.grid_S = S;
        result.grid_V = V;
        result.grid_payoff = payoff_vec;

        return result;
    }

    // ========================================================================
    // Barrier Option FD
    // ========================================================================

    FDResult fd_price_barrier(const BarrierOption &option,
                              const BlackScholesModel &model,
                              double S0,
                              const FDSettings &settings)
    {
        double K = option.get_strike();
        double T = option.get_maturity();
        bool is_call = option.is_call();
        double B = option.get_barrier();
        BarrierType type = option.get_barrier_type();

        // For knock-in options, use in-out parity: In = Vanilla - Out
        if (option.is_knock_in())
        {
            // Create corresponding vanilla option
            VanillaOption vanilla(K, T, is_call);
            FDResult vanilla_result = fd_price_european(vanilla, model, S0, settings);

            // Create corresponding knock-out barrier
            BarrierType out_type = (type == BarrierType::UpAndIn) ? BarrierType::UpAndOut : BarrierType::DownAndOut;
            BarrierOption out_option(K, T, is_call, out_type, B, option.get_rebate());
            FDResult out_result = fd_price_barrier(out_option, model, S0, settings);

            // In = Vanilla - Out
            FDResult result;
            result.price = vanilla_result.price - out_result.price;
            result.grid_S = vanilla_result.grid_S;
            result.grid_V.resize(vanilla_result.grid_V.size());
            for (size_t i = 0; i < vanilla_result.grid_V.size(); ++i)
            {
                result.grid_V[i] = vanilla_result.grid_V[i] - out_result.grid_V[i];
            }
            result.grid_payoff = vanilla_result.grid_payoff;
            return result;
        }

        // Knock-out option pricing with absorbing barrier
        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        double S_max = settings.get_S_max();
        int N_full = settings.get_N();
        int M = settings.get_M();
        double theta = settings.get_theta();

        // Determine grid domain based on barrier type
        double S_min_domain, S_max_domain;
        if (type == BarrierType::UpAndOut)
        {
            S_min_domain = 0.0;
            S_max_domain = B; // Domain [0, B]
        }
        else
        { // DownAndOut
            S_min_domain = B;
            S_max_domain = S_max; // Domain [B, S_max]
        }

        // Check if S0 is in valid domain
        if (type == BarrierType::UpAndOut && S0 >= B)
        {
            // Already knocked out
            FDResult result;
            result.price = option.get_rebate();
            result.grid_S = uniform_grid(S_max, N_full);
            result.grid_V.resize(N_full, 0.0);
            result.grid_payoff.resize(N_full, 0.0);
            return result;
        }
        if (type == BarrierType::DownAndOut && S0 <= B)
        {
            // Already knocked out
            FDResult result;
            result.price = option.get_rebate();
            result.grid_S = uniform_grid(S_max, N_full);
            result.grid_V.resize(N_full, 0.0);
            result.grid_payoff.resize(N_full, 0.0);
            return result;
        }

        // Build grid for the barrier domain
        int N = N_full;
        std::vector<double> S(N);
        double dS = (S_max_domain - S_min_domain) / (N - 1);
        for (int i = 0; i < N; ++i)
        {
            S[i] = S_min_domain + i * dS;
        }
        double dt = T / M;

        // Initialize option values at maturity
        std::vector<double> V(N);
        std::vector<double> payoff_vec(N);
        for (int i = 0; i < N; ++i)
        {
            V[i] = option.payoff(S[i]);
            payoff_vec[i] = V[i];
        }

        // At barrier, value = rebate (typically 0)
        if (type == BarrierType::UpAndOut)
        {
            V[N - 1] = option.get_rebate(); // S = B (upper boundary)
        }
        else
        {                               // DownAndOut
            V[0] = option.get_rebate(); // S = B (lower boundary)
        }

        // Build coefficient matrices
        std::vector<double> a_impl, b_impl, c_impl;
        std::vector<double> a_expl, b_expl, c_expl;
        build_tridiagonal_coeffs(S, r, q, sigma, dt, theta,
                                 a_impl, b_impl, c_impl,
                                 a_expl, b_expl, c_expl);

        // Time stepping
        std::vector<double> rhs(N);
        std::vector<double> V_new(N);

        for (int m = M - 1; m >= 0; --m)
        {
            double t = m * dt;

            // Apply explicit operator
            apply_explicit(a_expl, b_expl, c_expl, V, rhs);

            // Apply boundary conditions
            if (type == BarrierType::UpAndOut)
            {
                // S = 0: Call -> 0, Put -> K*exp(-r*t)
                if (is_call)
                {
                    rhs[0] = 0.0;
                }
                else
                {
                    rhs[0] = K * std::exp(-r * t);
                }
                // S = B: absorbing boundary (knocked out)
                rhs[N - 1] = option.get_rebate();
            }
            else
            { // DownAndOut
                // S = B: absorbing boundary
                rhs[0] = option.get_rebate();
                // S = S_max: Call -> S_max*exp(-q*t) - K*exp(-r*t), Put -> 0
                if (is_call)
                {
                    rhs[N - 1] = std::max(S_max_domain * std::exp(-q * t) - K * std::exp(-r * t), 0.0);
                }
                else
                {
                    rhs[N - 1] = 0.0;
                }
            }

            // Solve implicit system
            solve_tridiagonal(a_impl, b_impl, c_impl, rhs, V_new);

            // Enforce boundary conditions
            if (type == BarrierType::UpAndOut)
            {
                if (is_call)
                {
                    V_new[0] = 0.0;
                }
                else
                {
                    V_new[0] = K * std::exp(-r * t);
                }
                V_new[N - 1] = option.get_rebate();
            }
            else
            {
                V_new[0] = option.get_rebate();
                if (is_call)
                {
                    V_new[N - 1] = std::max(S_max_domain * std::exp(-q * t) - K * std::exp(-r * t), 0.0);
                }
                else
                {
                    V_new[N - 1] = 0.0;
                }
            }

            V = V_new;
        }

        // Interpolate to get V(S0)
        double price = linear_interp(S, V, S0);

        // Build full grid for output (with zeros outside domain)
        FDResult result;
        result.price = price;
        result.grid_S = uniform_grid(S_max, N_full);
        result.grid_V.resize(N_full, 0.0);
        result.grid_payoff.resize(N_full, 0.0);

        // Fill in values from computed domain
        for (int i = 0; i < N_full; ++i)
        {
            double s_val = result.grid_S[i];
            if (s_val >= S_min_domain && s_val <= S_max_domain)
            {
                result.grid_V[i] = linear_interp(S, V, s_val);
                result.grid_payoff[i] = option.payoff(s_val);
            }
        }

        return result;
    }

    // ========================================================================
    // Convenience functions
    // ========================================================================

    double fd_price(const VanillaOption &option,
                    const BlackScholesModel &model,
                    double S0,
                    const FDSettings &settings)
    {
        return fd_price_european(option, model, S0, settings).price;
    }

    double fd_price(const AmericanOption &option,
                    const BlackScholesModel &model,
                    double S0,
                    const FDSettings &settings)
    {
        return fd_price_american(option, model, S0, settings).price;
    }

    double fd_price(const BarrierOption &option,
                    const BlackScholesModel &model,
                    double S0,
                    const FDSettings &settings)
    {
        return fd_price_barrier(option, model, S0, settings).price;
    }

} // namespace OptionPricer
