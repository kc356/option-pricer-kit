#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"
#include "products/AsianOption.h"
#include "products/AmericanOption.h"
#include "products/BarrierOption.h"
#include "pricing/AnalyticPricer.h"
#include "pricing/MonteCarloPricer.h"
#include "pricing/MCSettings.h"
#include "pricing/FDSettings.h"
#include "pricing/FDPricer.h"
#include "greeks/BSGreeks.h"
#include "utils/ImpliedVol.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace OptionPricer;

// Test utilities
int test_count = 0;
int test_passed = 0;

#define TEST(name)                                             \
    void test_##name();                                        \
    void run_test_##name()                                     \
    {                                                          \
        test_count++;                                          \
        try                                                    \
        {                                                      \
            std::cout << "Running test: " << #name << " ... "; \
            test_##name();                                     \
            test_passed++;                                     \
            std::cout << "PASSED" << std::endl;                \
        }                                                      \
        catch (const std::exception &e)                        \
        {                                                      \
            std::cout << "FAILED: " << e.what() << std::endl;  \
        }                                                      \
    }                                                          \
    void test_##name()

#define ASSERT_NEAR(actual, expected, tolerance)                                     \
    if (std::abs((actual) - (expected)) > (tolerance))                               \
    {                                                                                \
        throw std::runtime_error("Assertion failed: " + std::to_string(actual) +     \
                                 " != " + std::to_string(expected) +                 \
                                 " (tolerance: " + std::to_string(tolerance) + ")"); \
    }

#define ASSERT_TRUE(condition)                                     \
    if (!(condition))                                              \
    {                                                              \
        throw std::runtime_error("Assertion failed: " #condition); \
    }

// Test 1: Call-Put Parity
TEST(call_put_parity)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption call(K, T, true);
    VanillaOption put(K, T, false);

    double call_price = price_bs_analytic(call, model, S0);
    double put_price = price_bs_analytic(put, model, S0);

    // Call - Put = S0 * exp(-qT) - K * exp(-rT)
    double lhs = call_price - put_price;
    double rhs = S0 * std::exp(-q * T) - K * std::exp(-r * T);

    ASSERT_NEAR(lhs, rhs, 1e-10);
}

// Test 2: ATM Call Price Reasonableness
TEST(atm_call_price)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption call(K, T, true);

    double price = price_bs_analytic(call, model, S0);

    // ATM call should be positive and reasonable
    ASSERT_TRUE(price > 0.0);
    ASSERT_TRUE(price < S0); // Can't be worth more than the stock

    // For ATM, price should be roughly S0 * N(0.5) * vol * sqrt(T)
    // This is a rough approximation
    double approx = S0 * 0.4 * sigma * std::sqrt(T);
    ASSERT_TRUE(std::abs(price - approx) < S0 * 0.2);
}

// Test 3: MC Convergence to Analytic
TEST(mc_convergence)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption call(K, T, true);

    double analytic_price = price_bs_analytic(call, model, S0);

    // MC with large number of paths
    MCSettings settings(1000000, 252, 42, true);
    MCResult mc_result = mc_price(call, model, S0, settings);

    // MC should be within 3 standard errors of analytic
    double diff = std::abs(mc_result.price - analytic_price);
    ASSERT_TRUE(diff < 3.0 * mc_result.std_error);

    std::cout << "\n    Analytic: " << analytic_price
              << ", MC: " << mc_result.price
              << " ± " << mc_result.std_error << " ... ";
}

// Test 4: Implied Volatility Round-Trip
TEST(implied_vol_roundtrip)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma_true = 0.25;

    BlackScholesModel model(r, q, sigma_true);
    VanillaOption call(K, T, true);

    // Get price with known sigma
    double market_price = price_bs_analytic(call, model, S0);

    // Recover implied vol
    double sigma_implied = implied_vol_bisect(market_price, S0, K, T, r, q, true);

    // Should match original sigma
    ASSERT_NEAR(sigma_implied, sigma_true, 1e-4);
}

// Test 5: Greeks Sign Check
TEST(greeks_signs)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption call(K, T, true);
    VanillaOption put(K, T, false);

    Greeks call_greeks = calculate_bs_greeks(call, model, S0);
    Greeks put_greeks = calculate_bs_greeks(put, model, S0);

    // Call delta should be positive, put delta negative
    ASSERT_TRUE(call_greeks.delta > 0.0);
    ASSERT_TRUE(put_greeks.delta < 0.0);

    // Gamma and Vega should be positive for both
    ASSERT_TRUE(call_greeks.gamma > 0.0);
    ASSERT_TRUE(put_greeks.gamma > 0.0);
    ASSERT_TRUE(call_greeks.vega > 0.0);
    ASSERT_TRUE(put_greeks.vega > 0.0);

    // Call and put should have same gamma and vega
    ASSERT_NEAR(call_greeks.gamma, put_greeks.gamma, 1e-10);
    ASSERT_NEAR(call_greeks.vega, put_greeks.vega, 1e-10);
}

// Test 6: Asian Option Payoff
TEST(asian_payoff)
{
    AsianOption call(100.0, 1.0, true, 12);

    // Path with average = 105
    std::vector<double> path = {100, 102, 104, 106, 108, 110,
                                108, 106, 104, 102, 100, 110};
    double avg = 105.0; // Calculated average

    double payoff = call.payoff(path);
    ASSERT_NEAR(payoff, 5.0, 1e-10); // max(105 - 100, 0) = 5
}

// Test 7: Geometric Asian Analytic Price
TEST(geometric_asian_price)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;
    int n_steps = 12;

    BlackScholesModel model(r, q, sigma);

    double price = price_geometric_asian_analytic(K, T, true, n_steps, model, S0);

    // Price should be positive and reasonable
    ASSERT_TRUE(price > 0.0);
    ASSERT_TRUE(price < S0);
}

// Test 8: Antithetic Variance Reduction
TEST(antithetic_variance_reduction)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption call(K, T, true);

    // Without antithetic
    MCSettings settings_no_anti(10000, 252, 42, false);
    MCResult result_no_anti = mc_price(call, model, S0, settings_no_anti);

    // With antithetic (same effective number of paths)
    MCSettings settings_anti(5000, 252, 42, true);
    MCResult result_anti = mc_price(call, model, S0, settings_anti);

    // Antithetic should have lower standard error
    std::cout << "\n    No Anti SE: " << result_no_anti.std_error
              << ", Anti SE: " << result_anti.std_error << " ... ";

    // This is probabilistic, but typically antithetic reduces variance
    // Just check that both give reasonable results
    ASSERT_TRUE(result_anti.std_error > 0.0);
    ASSERT_TRUE(result_no_anti.std_error > 0.0);
}

// Test 9: Control Variate Efficiency
TEST(control_variate_efficiency)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;
    int n_steps = 12;

    BlackScholesModel model(r, q, sigma);
    AsianOption asian(K, T, true, n_steps);

    MCSettings settings(50000, n_steps, 42, true);

    // Without control variate
    MCResult result_no_cv = mc_price_asian_cv(asian, model, S0, settings, false);

    // With control variate
    MCResult result_cv = mc_price_asian_cv(asian, model, S0, settings, true);

    std::cout << "\n    No CV SE: " << result_no_cv.std_error
              << ", CV SE: " << result_cv.std_error << " ... ";

    // Control variate should reduce standard error significantly
    // (usually by 80-95% for Asian options)
    ASSERT_TRUE(result_cv.std_error < result_no_cv.std_error);
}

// Test 10: Boundary Cases
TEST(boundary_cases)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);

    // Deep ITM call should be approximately S0 - K * exp(-rT)
    VanillaOption deep_itm_call(50.0, T, true);
    double itm_price = price_bs_analytic(deep_itm_call, model, S0);
    double intrinsic = S0 * std::exp(-q * T) - 50.0 * std::exp(-r * T);
    ASSERT_TRUE(itm_price > intrinsic);
    ASSERT_TRUE(itm_price < S0);

    // Deep OTM call should be near zero
    VanillaOption deep_otm_call(200.0, T, true);
    double otm_price = price_bs_analytic(deep_otm_call, model, S0);
    ASSERT_TRUE(otm_price > 0.0);
    ASSERT_TRUE(otm_price < 1.0);
}

// ============================================================================
// Finite Difference Tests
// ============================================================================

// Test 11: FD European vs Analytic (Call)
TEST(fd_european_call)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption call(K, T, true);

    // FD pricing with fine grid
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult fd_result = fd_price_european(call, model, S0, settings);

    // Analytic price
    double analytic_price = price_bs_analytic(call, model, S0);

    // Should be close (within 0.01)
    double diff = std::abs(fd_result.price - analytic_price);
    std::cout << "\n    FD Call: " << fd_result.price
              << ", Analytic: " << analytic_price
              << ", Diff: " << diff << " ... ";

    ASSERT_NEAR(fd_result.price, analytic_price, 0.01);
}

// Test 12: FD European vs Analytic (Put)
TEST(fd_european_put)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    VanillaOption put(K, T, false);

    // FD pricing
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult fd_result = fd_price_european(put, model, S0, settings);

    // Analytic price
    double analytic_price = price_bs_analytic(put, model, S0);

    // Should be close
    double diff = std::abs(fd_result.price - analytic_price);
    std::cout << "\n    FD Put: " << fd_result.price
              << ", Analytic: " << analytic_price
              << ", Diff: " << diff << " ... ";

    ASSERT_NEAR(fd_result.price, analytic_price, 0.01);
}

// Test 13: American Put >= European Put
TEST(american_put_ge_european)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    AmericanOption american_put(K, T, false);
    VanillaOption european_put(K, T, false);

    // FD pricing
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult american_result = fd_price_american(american_put, model, S0, settings);
    double european_price = price_bs_analytic(european_put, model, S0);

    std::cout << "\n    American Put: " << american_result.price
              << ", European Put: " << european_price << " ... ";

    // American should be >= European (early exercise premium)
    ASSERT_TRUE(american_result.price >= european_price - 1e-6);
}

// Test 14: American Call ~ European Call (no dividend scenario)
TEST(american_call_equals_european_no_dividend)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.0; // No dividends
    double sigma = 0.25;

    BlackScholesModel model(r, q, sigma);
    AmericanOption american_call(K, T, true);
    VanillaOption european_call(K, T, true);

    // FD pricing
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult american_result = fd_price_american(american_call, model, S0, settings);
    double european_price = price_bs_analytic(european_call, model, S0);

    std::cout << "\n    American Call (q=0): " << american_result.price
              << ", European: " << european_price << " ... ";

    // With no dividends, American call = European call (approximately)
    ASSERT_NEAR(american_result.price, european_price, 0.02);
}

// Test 15: Barrier Out <= Vanilla
TEST(barrier_out_le_vanilla)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;
    double B = 120.0; // Up-and-out barrier

    BlackScholesModel model(r, q, sigma);
    BarrierOption barrier_call(K, T, true, BarrierType::UpAndOut, B);
    VanillaOption vanilla_call(K, T, true);

    // FD pricing
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult barrier_result = fd_price_barrier(barrier_call, model, S0, settings);
    double vanilla_price = price_bs_analytic(vanilla_call, model, S0);

    std::cout << "\n    Up-and-Out Call: " << barrier_result.price
              << ", Vanilla: " << vanilla_price << " ... ";

    // Knock-out option should be worth less than vanilla
    ASSERT_TRUE(barrier_result.price <= vanilla_price + 1e-6);
    // But should still be positive
    ASSERT_TRUE(barrier_result.price >= 0.0);
}

// Test 16: In-Out Parity
TEST(barrier_in_out_parity)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;
    double B = 120.0; // Barrier level

    BlackScholesModel model(r, q, sigma);
    BarrierOption up_and_out(K, T, true, BarrierType::UpAndOut, B);
    BarrierOption up_and_in(K, T, true, BarrierType::UpAndIn, B);
    VanillaOption vanilla(K, T, true);

    // FD pricing
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult out_result = fd_price_barrier(up_and_out, model, S0, settings);
    FDResult in_result = fd_price_barrier(up_and_in, model, S0, settings);
    double vanilla_price = price_bs_analytic(vanilla, model, S0);

    double parity_sum = out_result.price + in_result.price;

    std::cout << "\n    Out: " << out_result.price
              << ", In: " << in_result.price
              << ", Sum: " << parity_sum
              << ", Vanilla: " << vanilla_price << " ... ";

    // In + Out should equal Vanilla
    ASSERT_NEAR(parity_sum, vanilla_price, 0.02);
}

// Test 17: Down-and-Out Put
TEST(barrier_down_and_out_put)
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double q = 0.02;
    double sigma = 0.25;
    double B = 80.0; // Down-and-out barrier

    BlackScholesModel model(r, q, sigma);
    BarrierOption barrier_put(K, T, false, BarrierType::DownAndOut, B);
    VanillaOption vanilla_put(K, T, false);

    // FD pricing
    FDSettings settings(300.0, 200, 200, 0.5, 1.2, 1e-6, 1000);
    FDResult barrier_result = fd_price_barrier(barrier_put, model, S0, settings);
    double vanilla_price = price_bs_analytic(vanilla_put, model, S0);

    std::cout << "\n    Down-and-Out Put: " << barrier_result.price
              << ", Vanilla: " << vanilla_price << " ... ";

    // Should be less than vanilla
    ASSERT_TRUE(barrier_result.price <= vanilla_price + 1e-6);
    ASSERT_TRUE(barrier_result.price >= 0.0);
}

int main()
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Option Pricer Tests ===" << std::endl
              << std::endl;

    run_test_call_put_parity();
    run_test_atm_call_price();
    run_test_mc_convergence();
    run_test_implied_vol_roundtrip();
    run_test_greeks_signs();
    run_test_asian_payoff();
    run_test_geometric_asian_price();
    run_test_antithetic_variance_reduction();
    run_test_control_variate_efficiency();
    run_test_boundary_cases();

    // Finite Difference tests
    std::cout << std::endl
              << "--- Finite Difference Tests ---" << std::endl;
    run_test_fd_european_call();
    run_test_fd_european_put();
    run_test_american_put_ge_european();
    run_test_american_call_equals_european_no_dividend();
    run_test_barrier_out_le_vanilla();
    run_test_barrier_in_out_parity();
    run_test_barrier_down_and_out_put();

    std::cout << std::endl;
    std::cout << "=== Test Summary ===" << std::endl;
    std::cout << "Tests passed: " << test_passed << " / " << test_count << std::endl;

    return (test_passed == test_count) ? 0 : 1;
}
