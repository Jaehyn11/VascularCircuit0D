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

TEST_CASE("HANDLE INFLOW CSV", "[input, flows]") {

  InflowData inflow_data;

  SECTION("Non-existent file") {
    std::string inflow_csv = "non_existent_file.csv";
    bool result = load_inflow_csv(inflow_csv, inflow_data);
    REQUIRE(result == false);
  }
  SECTION("Existing file") {
    std::string inflow_csv = "input/AAo.csv";
    bool result = load_inflow_csv(inflow_csv, inflow_data);
    REQUIRE(result == true);
    REQUIRE(inflow_data.loaded == true);
    REQUIRE(inflow_data.t.size() == 40);
    REQUIRE(inflow_data.q.size() == 40);

    // Checking the 14th point (index 13) because I don't want to hardcode the
    // whole file
    REQUIRE_THAT(inflow_data.t[13],
                 Catch::Matchers::WithinAbs(0.338755980861244, 1e-6));
    REQUIRE_THAT(inflow_data.q[13],
                 Catch::Matchers::WithinAbs(18.04337363041816, 1e-6));
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

TEST_CASE("Q_IN_FROM_CSV", "[input, flows]") {
  InflowData inflow_data;
  std::string inflow_csv = "input/AAo.csv";
  bool result = load_inflow_csv(inflow_csv, inflow_data);

  SECTION("Interpolation") {
    double t = 0.0;
    double q_in = Q_in_from_csv(t, inflow_data);
    REQUIRE_THAT(q_in, Catch::Matchers::WithinAbs(-2.67291413040408798, 1e-6));
  }

  SECTION("Clamping") {
    double t = 1.0133971291866029;
    double q_in = Q_in_from_csv(t, inflow_data);
    REQUIRE_THAT(q_in, Catch::Matchers::WithinAbs(-4.733727810650976, 1e-6));
  }
}

TEST_CASE("Q_IN", "[flows]") {
  Parameters par;
  par.T_cycle = 1.0;
  double t = 0.0;

  InflowData inflow_data;
  std::string inflow_csv = "input/AAo.csv";
  bool result = load_inflow_csv(inflow_csv, inflow_data);

  SECTION("Analytic") {
    bool use_csv = false;
    double q_in = Q_in(t, par, inflow_data, use_csv);
    REQUIRE_THAT(q_in, Catch::Matchers::WithinAbs(80.0, 1e-6));
  }

  SECTION("CSV") {
    bool use_csv = true;
    double q_in = Q_in(t, par, inflow_data, use_csv);
    REQUIRE_THAT(q_in, Catch::Matchers::WithinAbs(-2.67291413040408798, 1e-6));
  }
}

TEST_CASE("SOLVER CHECKS", "[solver]") {
  double t = 0.0;
  Parameters par;
  State s0;
  State s1;
  double t_final, dt;
  bool use_csv;
  std::string inflow_csv;

  bool result1 = load_input_file("input/input.inp", par, s0, t_final, dt,
                                 use_csv, inflow_csv);

  InflowData inflow_data;
  bool result2 = load_inflow_csv(inflow_csv, inflow_data);

  rhs_full(t, s0, par, inflow_data, use_csv, s1);

  REQUIRE_THAT(s1.Q[0], Catch::Matchers::WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(s1.Q[4], Catch::Matchers::WithinAbs(0.0, 1e-6));
}