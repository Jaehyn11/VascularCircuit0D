#include "solver_core.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("BASIC SANITY TEST", "[sanity]") { REQUIRE(1 + 1 == 2); }

// This is a single test case but with multiple parts
// Each SECTION is not counted as a separate test case
// in the command line (for reference)
TEST_CASE("HANDLE INPUT FILE CORRECTLY", "[input]") {
  Parameters par;
  State s0;
  double t_final, dt;
  bool use_csv;
  std::string inflow_csv;

  SECTION("Non-existent file") {
    bool result = load_input_file("non_existent_file.inp", par, s0, t_final, dt,
                                  use_csv, inflow_csv);
    REQUIRE(result == false);
  }
  // NOTE THAT THIS IS BASED ON THE INPUT FILE COPIED DURING BUILD...
  // THAT IS, THE TEST WILL FAIL IF THE COPIED INPUT FILE IS CHANGED AFTER
  // BUILDING
  SECTION("Existing file") {
    constexpr int n_nodes = 5;
    constexpr int n_segs = 5;

    std::array<double, n_nodes> C_expected{1.47, 0.49, 0.088, 0.0147, 0.0};
    std::array<double, n_segs> R_expected{0.15, 0.15, 0.35, 0.35, 0.00};
    std::array<double, n_segs> L_expected{0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<double, n_nodes> P_initial_expected{75.0, 75.0, 75.0, 75.0, 0.0};

    bool result = load_input_file("input/input.inp", par, s0, t_final, dt,
                                  use_csv, inflow_csv);
    // Load is treated as successful
    REQUIRE(result == true);
    // Check if the parameters are loaded correctly
    // Unfortunately, Catch2 does not appear to have a built-in for std::array
    // although it does have one for std::vector... damn (so need to loop)
    for (size_t i = 0; i < n_nodes; ++i) {
      REQUIRE_THAT(par.C[i], Catch::Matchers::WithinAbs(C_expected[i], 1e-6));
    }
    for (size_t i = 0; i < n_segs; ++i) {
      REQUIRE_THAT(par.R[i], Catch::Matchers::WithinAbs(R_expected[i], 1e-6));
      REQUIRE_THAT(par.L[i], Catch::Matchers::WithinAbs(L_expected[i], 1e-6));
    }
    REQUIRE_THAT(par.P_RA, Catch::Matchers::WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(par.T_cycle, Catch::Matchers::WithinAbs(1.0, 1e-6));

    // Check initial state
    for (size_t i = 0; i < n_nodes; ++i) {
      REQUIRE_THAT(s0.P[i],
                   Catch::Matchers::WithinAbs(P_initial_expected[i], 1e-6));
    }
    // Flows are initialized with separate function... maybe check that
    // separately first?
  }
}

TEST_CASE("ANALYTIC FLOWS", "[flows]") {
  Parameters par;
  par.T_cycle = 1.0;

  SECTION("Zero") {
    double t = 0.0;
    double flow = Q_in_analytic(t, par);
    REQUIRE_THAT(flow, Catch::Matchers::WithinAbs(80.0, 1e-6));
  }
  SECTION("Peak") {
    double t = 0.25;
    double flow = Q_in_analytic(t, par);
    REQUIRE_THAT(flow, Catch::Matchers::WithinAbs(120.0, 1e-6));
  }
}