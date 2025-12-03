#include "solver_core.hpp"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("BASIC SANITY TEST", "[sanity]") { REQUIRE(1 + 1 == 2); }

TEST_CASE("HANDLE INPUT FILE ERROR", "[input, error]") {
  Parameters par;
  State s0;
  double t_final, dt;
  bool use_csv;
  std::string inflow_csv;
  bool result = load_input_file("non_existent_file.inp", par, s0, t_final, dt,
                                use_csv, inflow_csv);
  REQUIRE(result == false);
}