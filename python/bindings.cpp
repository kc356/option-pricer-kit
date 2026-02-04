#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

namespace py = pybind11;
using namespace OptionPricer;

PYBIND11_MODULE(option_pricer_cpp, m)
{
     m.doc() = "Option pricing library with C++ core";

     // BlackScholesModel
     py::class_<BlackScholesModel>(m, "BlackScholesModel")
         .def(py::init<double, double, double>(),
              py::arg("r"), py::arg("q"), py::arg("sigma"),
              "Black-Scholes model with constant volatility")
         .def_property("r", &BlackScholesModel::get_r, &BlackScholesModel::set_r,
                       "Risk-free rate")
         .def_property("q", &BlackScholesModel::get_q, &BlackScholesModel::set_q,
                       "Dividend yield")
         .def_property("sigma", &BlackScholesModel::get_sigma, &BlackScholesModel::set_sigma,
                       "Volatility")
         .def("__repr__", [](const BlackScholesModel &m)
              { return "<BlackScholesModel r=" + std::to_string(m.get_r()) +
                       " q=" + std::to_string(m.get_q()) +
                       " sigma=" + std::to_string(m.get_sigma()) + ">"; });

     // VanillaOption
     py::class_<VanillaOption>(m, "VanillaOption")
         .def(py::init<double, double, bool>(),
              py::arg("K"), py::arg("T"), py::arg("is_call"),
              "European vanilla option")
         .def_property_readonly("K", &VanillaOption::get_strike, "Strike price")
         .def_property_readonly("T", &VanillaOption::get_maturity, "Time to maturity")
         .def_property_readonly("is_call", &VanillaOption::is_call, "Call or put")
         .def("payoff", &VanillaOption::payoff, py::arg("S_T"),
              "Calculate payoff at maturity")
         .def("__repr__", [](const VanillaOption &o)
              { return "<VanillaOption " + std::string(o.is_call() ? "Call" : "Put") +
                       " K=" + std::to_string(o.get_strike()) +
                       " T=" + std::to_string(o.get_maturity()) + ">"; });

     // AsianOption
     py::class_<AsianOption>(m, "AsianOption")
         .def(py::init<double, double, bool, int>(),
              py::arg("K"), py::arg("T"), py::arg("is_call"), py::arg("n_steps"),
              "Asian option with arithmetic average")
         .def_property_readonly("K", &AsianOption::get_strike, "Strike price")
         .def_property_readonly("T", &AsianOption::get_maturity, "Time to maturity")
         .def_property_readonly("is_call", &AsianOption::is_call, "Call or put")
         .def_property_readonly("n_steps", &AsianOption::get_n_steps, "Number of monitoring points")
         .def("payoff", &AsianOption::payoff, py::arg("path"),
              "Calculate payoff from arithmetic average")
         .def("payoff_geometric", &AsianOption::payoff_geometric, py::arg("path"),
              "Calculate payoff from geometric average")
         .def("__repr__", [](const AsianOption &o)
              { return "<AsianOption " + std::string(o.is_call() ? "Call" : "Put") +
                       " K=" + std::to_string(o.get_strike()) +
                       " T=" + std::to_string(o.get_maturity()) +
                       " n_steps=" + std::to_string(o.get_n_steps()) + ">"; });

     // AmericanOption
     py::class_<AmericanOption>(m, "AmericanOption")
         .def(py::init<double, double, bool>(),
              py::arg("K"), py::arg("T"), py::arg("is_call"),
              "American vanilla option with early exercise")
         .def_property_readonly("K", &AmericanOption::get_strike, "Strike price")
         .def_property_readonly("T", &AmericanOption::get_maturity, "Time to maturity")
         .def_property_readonly("is_call", &AmericanOption::is_call, "Call or put")
         .def("payoff", &AmericanOption::payoff, py::arg("S"),
              "Calculate payoff at spot price S")
         .def("exercise_value", &AmericanOption::exercise_value, py::arg("S"),
              "Calculate immediate exercise value")
         .def("__repr__", [](const AmericanOption &o)
              { return "<AmericanOption " + std::string(o.is_call() ? "Call" : "Put") +
                       " K=" + std::to_string(o.get_strike()) +
                       " T=" + std::to_string(o.get_maturity()) + ">"; });

     // BarrierType enum
     py::enum_<BarrierType>(m, "BarrierType")
         .value("UpAndOut", BarrierType::UpAndOut)
         .value("DownAndOut", BarrierType::DownAndOut)
         .value("UpAndIn", BarrierType::UpAndIn)
         .value("DownAndIn", BarrierType::DownAndIn)
         .export_values();

     // BarrierOption
     py::class_<BarrierOption>(m, "BarrierOption")
         .def(py::init<double, double, bool, BarrierType, double, double>(),
              py::arg("K"), py::arg("T"), py::arg("is_call"),
              py::arg("barrier_type"), py::arg("B"),
              py::arg("rebate") = 0.0,
              "European barrier option")
         .def_property_readonly("K", &BarrierOption::get_strike, "Strike price")
         .def_property_readonly("T", &BarrierOption::get_maturity, "Time to maturity")
         .def_property_readonly("is_call", &BarrierOption::is_call, "Call or put")
         .def_property_readonly("barrier_type", &BarrierOption::get_barrier_type, "Barrier type")
         .def_property_readonly("B", &BarrierOption::get_barrier, "Barrier level")
         .def_property_readonly("rebate", &BarrierOption::get_rebate, "Rebate on knock-out")
         .def("is_knock_out", &BarrierOption::is_knock_out, "Check if knock-out type")
         .def("is_knock_in", &BarrierOption::is_knock_in, "Check if knock-in type")
         .def("is_up", &BarrierOption::is_up, "Check if up-type barrier")
         .def("is_down", &BarrierOption::is_down, "Check if down-type barrier")
         .def("payoff", &BarrierOption::payoff, py::arg("S_T"),
              "Calculate payoff at maturity (if not knocked out)")
         .def("is_breached", &BarrierOption::is_breached, py::arg("S"),
              "Check if barrier has been breached")
         .def("__repr__", [](const BarrierOption &o)
              { return "<BarrierOption " + std::string(o.is_call() ? "Call" : "Put") +
                       " K=" + std::to_string(o.get_strike()) +
                       " T=" + std::to_string(o.get_maturity()) +
                       " B=" + std::to_string(o.get_barrier()) +
                       " " + std::string(barrier_type_to_string(o.get_barrier_type())) + ">"; });

     // MCSettings
     py::class_<MCSettings>(m, "MCSettings")
         .def(py::init<int, int, unsigned int, bool>(),
              py::arg("n_paths") = 100000,
              py::arg("n_steps") = 252,
              py::arg("seed") = 42,
              py::arg("antithetic") = true,
              "Monte Carlo simulation settings")
         .def_property("n_paths", &MCSettings::get_n_paths, &MCSettings::set_n_paths,
                       "Number of simulation paths")
         .def_property("n_steps", &MCSettings::get_n_steps, &MCSettings::set_n_steps,
                       "Number of time steps")
         .def_property("seed", &MCSettings::get_seed, &MCSettings::set_seed,
                       "Random seed")
         .def_property("antithetic", &MCSettings::use_antithetic, &MCSettings::set_antithetic,
                       "Use antithetic variates")
         .def("__repr__", [](const MCSettings &s)
              { return "<MCSettings n_paths=" + std::to_string(s.get_n_paths()) +
                       " n_steps=" + std::to_string(s.get_n_steps()) +
                       " seed=" + std::to_string(s.get_seed()) +
                       " antithetic=" + (s.use_antithetic() ? "True" : "False") + ">"; });

     // FDSettings
     py::class_<FDSettings>(m, "FDSettings")
         .def(py::init<double, int, int, double, double, double, int>(),
              py::arg("S_max") = 300.0,
              py::arg("N") = 200,
              py::arg("M") = 200,
              py::arg("theta") = 0.5,
              py::arg("omega") = 1.2,
              py::arg("tol") = 1e-6,
              py::arg("max_iter") = 1000,
              "Finite Difference settings for PDE pricing")
         .def_property("S_max", &FDSettings::get_S_max, &FDSettings::set_S_max,
                       "Upper bound of S domain")
         .def_property("N", &FDSettings::get_N, &FDSettings::set_N,
                       "Number of space grid points")
         .def_property("M", &FDSettings::get_M, &FDSettings::set_M,
                       "Number of time steps")
         .def_property("theta", &FDSettings::get_theta, &FDSettings::set_theta,
                       "Scheme parameter (0.5 = Crank-Nicolson)")
         .def_property("omega", &FDSettings::get_omega, &FDSettings::set_omega,
                       "PSOR relaxation parameter")
         .def_property("tol", &FDSettings::get_tol, &FDSettings::set_tol,
                       "PSOR convergence tolerance")
         .def_property("max_iter", &FDSettings::get_max_iter, &FDSettings::set_max_iter,
                       "Maximum PSOR iterations")
         .def_property_readonly("dS", &FDSettings::get_dS,
                                "Grid spacing")
         .def("__repr__", [](const FDSettings &s)
              { return "<FDSettings S_max=" + std::to_string(s.get_S_max()) +
                       " N=" + std::to_string(s.get_N()) +
                       " M=" + std::to_string(s.get_M()) +
                       " theta=" + std::to_string(s.get_theta()) + ">"; });

     // FDResult
     py::class_<FDResult>(m, "FDResult")
         .def_readonly("price", &FDResult::price, "Option price at S0")
         .def_readonly("grid_S", &FDResult::grid_S, "Spot price grid")
         .def_readonly("grid_V", &FDResult::grid_V, "Option values at t=0")
         .def_readonly("grid_payoff", &FDResult::grid_payoff, "Payoff values")
         .def("__repr__", [](const FDResult &r)
              { return "<FDResult price=" + std::to_string(r.price) + ">"; });

     // MCResult
     py::class_<MCResult>(m, "MCResult")
         .def_readonly("price", &MCResult::price, "Estimated price")
         .def_readonly("std_error", &MCResult::std_error, "Standard error")
         .def_readonly("ci_lower", &MCResult::ci_lower, "95% CI lower bound")
         .def_readonly("ci_upper", &MCResult::ci_upper, "95% CI upper bound")
         .def("__repr__", [](const MCResult &r)
              { return "<MCResult price=" + std::to_string(r.price) +
                       " std_error=" + std::to_string(r.std_error) + ">"; });

     // Greeks
     py::class_<Greeks>(m, "Greeks")
         .def_readonly("delta", &Greeks::delta, "Delta: ∂V/∂S")
         .def_readonly("gamma", &Greeks::gamma, "Gamma: ∂²V/∂S²")
         .def_readonly("vega", &Greeks::vega, "Vega: ∂V/∂σ (per 1%)")
         .def_readonly("theta", &Greeks::theta, "Theta: ∂V/∂t")
         .def_readonly("rho", &Greeks::rho, "Rho: ∂V/∂r")
         .def("__repr__", [](const Greeks &g)
              { return "<Greeks delta=" + std::to_string(g.delta) +
                       " gamma=" + std::to_string(g.gamma) +
                       " vega=" + std::to_string(g.vega) + ">"; });

     // Pricing functions
     m.def("bs_price", &price_bs_analytic,
           py::arg("option"), py::arg("model"), py::arg("S0"),
           "Black-Scholes analytic price for European vanilla option");

     m.def("bs_greeks", &calculate_bs_greeks,
           py::arg("option"), py::arg("model"), py::arg("S0"),
           "Calculate Greeks for European vanilla option");

     // MC pricing - overloaded versions
     m.def("mc_price",
           py::overload_cast<const VanillaOption &, const BlackScholesModel &, double, const MCSettings &>(&mc_price),
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Monte Carlo price for vanilla option");

     m.def("mc_price",
           py::overload_cast<const AsianOption &, const BlackScholesModel &, double, const MCSettings &>(&mc_price),
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Monte Carlo price for Asian option");

     m.def("mc_price_asian_cv", &mc_price_asian_cv,
           py::arg("option"), py::arg("model"), py::arg("S0"),
           py::arg("settings"), py::arg("use_cv") = true,
           "Monte Carlo price for Asian option with control variate");

     m.def("geometric_asian_price", &price_geometric_asian_analytic,
           py::arg("K"), py::arg("T"), py::arg("is_call"),
           py::arg("n_steps"), py::arg("model"), py::arg("S0"),
           "Analytic price for geometric Asian option");

     m.def("implied_vol", &implied_vol_bisect,
           py::arg("market_price"), py::arg("S0"), py::arg("K"),
           py::arg("T"), py::arg("r"), py::arg("q"), py::arg("is_call"),
           py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
           "Calculate implied volatility using bisection method");

     // FD pricing functions
     m.def("fd_price_european", &fd_price_european,
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Finite Difference price for European vanilla option (returns FDResult)");

     m.def("fd_price_american", &fd_price_american,
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Finite Difference price for American vanilla option with PSOR (returns FDResult)");

     m.def("fd_price_barrier", &fd_price_barrier,
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Finite Difference price for European barrier option (returns FDResult)");

     // Overloaded fd_price convenience functions
     m.def("fd_price",
           py::overload_cast<const VanillaOption &, const BlackScholesModel &, double, const FDSettings &>(&fd_price),
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Finite Difference price for European vanilla option");

     m.def("fd_price",
           py::overload_cast<const AmericanOption &, const BlackScholesModel &, double, const FDSettings &>(&fd_price),
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Finite Difference price for American option");

     m.def("fd_price",
           py::overload_cast<const BarrierOption &, const BlackScholesModel &, double, const FDSettings &>(&fd_price),
           py::arg("option"), py::arg("model"), py::arg("S0"), py::arg("settings"),
           "Finite Difference price for barrier option");
}
