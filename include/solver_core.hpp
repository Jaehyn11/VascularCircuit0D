#pragma once

#include <array>
#include <string>
#include <vector>

constexpr int N_NODES = 5;
constexpr int N_SEGS = 5; // 0..3: between nodes, 4: last node -> RA

// Parameters and state used by solver core
// 1-1. Parameters for the simulation (theta vector)
struct Parameters {
  std::array<double, N_NODES> C; // Compliances [mL/mmHg] C[i] at node i
  std::array<double, N_SEGS> R;  // Segment resistances [mmHg·s/mL] for flows
  std::array<double, N_SEGS> L;  // Inductances: L[e] at segment e

  double P_RA;    // Right atrium pressure [mmHg]
  double T_cycle; // Cardiac period (for analytic inflow) [s]
};

// 1-2. State vector: node pressures (y vector)
struct State {
  std::array<double, N_NODES> P; // P[i] = pressure at node i
  std::array<double, N_SEGS> Q;  // flows in segments
};

// 2. Inflow data aquisition
struct InflowData {
  std::vector<double> t; // time [s]
  std::vector<double> q; // flow [mL/s]
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