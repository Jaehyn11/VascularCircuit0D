#include "solver_core.hpp"
#include <array>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>

//---------------------------------------------
// 0D Vascular circuit model: linear chain
// node 0 -> node 1 -> node 2 -> node 3 -> node 4 -> RA
// Units:
//   Q:     mL/s
//   P:     mmHg
//   C:     mL/mmHg
//   R:     mmHg·s/mL
//   t, dt: s
//---------------------------------------------

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " path/to/input.inp\n";
    return 1;
  }

  Parameters par;
  State s, k1, k2;
  double t_final, dt;
  bool use_csv;
  std::string input_file = argv[1];
  std::string inflow_csv_name;

  // 1. Read input.inp for theta, ICs, and time-stepping
  if (!load_input_file(input_file, par, s, t_final, dt, use_csv,
                       inflow_csv_name)) {
    return 1;
  }

  // 2. Inflow data from CSV
  InflowData inflow;
  if (use_csv) {
    // Build path from $HOME/NERSE570/VascularCircuit/input/<inflow_csv_name>
    // //setme
    const char *home_env = std::getenv("HOME");
    std::string inflow_file;
    if (home_env) {
      std::string home(home_env);
      // Use project-local directory under $HOME
      inflow_file = "./input/" + inflow_csv_name; // setme
    } else {
      inflow_file = inflow_csv_name;
    }

    std::cout << "Trying inflow file: " << inflow_file << "\n";
    if (!load_inflow_csv(inflow_file, inflow)) {
      std::cerr << "Falling back to analytic inflow.\n";
      use_csv = false;
    }
  }

  // 4. Time stepping setup
  int n_steps = static_cast<int>(t_final / dt);
  std::ofstream fout("output.csv");

  // output.csv header
  fout << "t";
  for (int i = 0; i < N_NODES; ++i)
    fout << ",P_" << i;
  for (int e = 0; e < N_SEGS; ++e)
    fout << ",Q_" << e;
  fout << "\n";

  double t = 0.0;

  // Compute Q for output (IMPORTANT: dynamic if L>0, Ohm's law if L==0)
  auto compute_Qout = [&](const State &s, const Parameters &par) {
    std::array<double, N_SEGS> Qout;
    for (size_t e = 0; e < N_SEGS; ++e) {
      double Pup, Pdown;
      if (e < N_NODES - 1) {
        Pup = s.P[e];
        Pdown = s.P[e + 1];
      } else {
        Pup = s.P[N_NODES - 1];
        Pdown = par.P_RA;
      }

      if (par.L[e] > 0.0) {
        Qout[e] = s.Q[e]; // dynamic segment
      } else {
        Qout[e] = (Pup - Pdown) / par.R[e]; // purely resistive
      }
    }
    return Qout;
  };

  // write initial state
  auto Qout0 = compute_Qout(s, par);
  fout << t;
  for (size_t i = 0; i < N_NODES; ++i)
    fout << "," << s.P[i];
  for (size_t e = 0; e < N_SEGS; ++e)
    fout << "," << Qout0[e];
  fout << "\n";

  // Diagnostics: flow through first and last segment
  double Q0_sum = 0.0, Qlast_sum = 0.0;
  int count = 0;

  // 5. RK2 loop
  // RK2 is chosen as a simple, 2nd-order accurate explicit scheme.
  // ------------------------------------------------------------
  // Time integrator: explicit RK2 (Heun method)
  // ------------------------------------------------------------
  // NOTE for future Julia Han:
  //   - This loop is the only place where the time integration
  //     scheme is hard-coded.
  //   - The physics of the 0D model is implemented in rhs_full(),
  //     which computes dP/dt and dQ/dt for a given state.
  //   - If you want to:
  //        * replace RK2 with RK4, BDF, etc., or
  //        * use PETSc TS with MPI / time-parallelization,
  //     you should replace this RK2 loop with your own time
  //     integrator, while reusing rhs_full() and the State/
  //     Parameters/InflowData structures.
  // ------------------------------------------------------------
  for (int n = 0; n < n_steps; ++n) {

    // Stage 1: evaluate RHS at current time/state (t, s)
    // k1 ~= dy/dt at the beginning of the step
    rhs_full(t, s, par, inflow, use_csv, k1);

    // Predictor step: forward Euler guess at t + dt
    // s_pred = s^n + dt * k1
    State s_pred;
    for (size_t i = 0; i < N_NODES; ++i) {
      s_pred.P[i] = s.P[i] + dt * k1.P[i]; // predicted pressures
    }
    for (size_t e = 0; e < N_SEGS; ++e) {
      s_pred.Q[e] = s.Q[e] + dt * k1.Q[e]; // predicted flows
    }

    // Stage 2: evaluate RHS at predicted state (t + dt, s_pred)
    // k2 ~= dy/dt at the end of the step
    rhs_full(t + dt, s_pred, par, inflow, use_csv, k2);

    // Corrector step: average slopes k1 and k2
    // s^{n+1} = s^n + (dt/2) * (k1 + k2)
    for (size_t i = 0; i < N_NODES; ++i) {
      s.P[i] += 0.5 * dt * (k1.P[i] + k2.P[i]);
    }
    for (size_t e = 0; e < N_SEGS; ++e) {
      s.Q[e] += 0.5 * dt * (k1.Q[e] + k2.Q[e]);
    }

    // advance physical time
    t += dt;

    // Compute flows for output/diagnostics
    // auto Qout = compute_Qout(s, par); // commented out because unused and
    // -Werror is on

    // output.csv header: time, pressures, flows for this step
    fout << t;
    for (size_t i = 0; i < N_NODES; ++i)
      fout << "," << s.P[i];
    for (size_t e = 0; e < N_SEGS; ++e)
      fout << "," << s.Q[e];
    fout << "\n";

    // Diagnostics: accumulate flows in segment 0 and last segment
    Q0_sum += s.Q[0];
    Qlast_sum += s.Q[N_SEGS - 1];
    count++;
  }

  fout.close();
  std::cout << "Simulation completed. Results written to output.csv\n";
  if (count > 0) {
    std::cout << "Mean Q_seg0  = " << Q0_sum / count << " mL/s\n";
    std::cout << "Mean Q_last  = " << Qlast_sum / count << " mL/s\n";
  }

  return 0;
}