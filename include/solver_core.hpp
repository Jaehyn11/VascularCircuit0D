#pragma once

#include <array>
#include <string>
#include <vector>

constexpr int N_NODES = 5;
constexpr int N_SEGS = 5;

// Parameters and state used by solver core
struct Parameters {
  std::array<double, N_NODES> C;
  std::array<double, N_SEGS> R;
  std::array<double, N_SEGS> L;
  double P_RA;
  double T_cycle;
};

struct State {
  std::array<double, N_NODES> P;
  std::array<double, N_SEGS> Q;
};

struct InflowData {
  std::vector<double> t;
  std::vector<double> q;
  bool loaded = false;
};

// I/O and helper functions used by tests and the main executable
bool load_input_file(const std::string &filename, Parameters &par, State &s0,
                     double &t_final, double &dt, bool &use_csv,
                     std::string &inflow_csv);

bool load_inflow_csv(const std::string &filename, InflowData &inflow);

double Q_in_analytic(double t, const Parameters &par);
double Q_in_from_csv(double t, const InflowData &inflow);
double Q_in(double t, const Parameters &par, const InflowData &inflow,
            bool use_csv);

void rhs_full(double t, const State &s, const Parameters &par,
              const InflowData &inflow, bool use_csv, State &ddt);
#pragma once
// ADD PROTOTYPES AS NEEDED FOR TESTS.