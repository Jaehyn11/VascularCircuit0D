#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

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

constexpr int N_NODES = 5;
constexpr int N_SEGS = 5; // 0..3: between nodes, 4: last node -> RA

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

// 1-3 Helper function to load VascularCircuit0D/input/input.inp file
bool load_input_file(const std::string &filename, Parameters &par, State &s0,
                     double &t_final, double &dt, bool &use_csv,
                     std::string &inflow_csv) {
  std::ifstream fin(filename);
  if (!fin) {
    std::cerr << "Error: cannot open input file: " << filename << "\n";
    return false;
  }

  // Default parameters
  par.C.fill(1.0); // mL/mmHg
  par.R.fill(1.0); // mmHg·s/mL
  par.L.fill(0.0); // inertia disabled by default (stable)
  par.P_RA = 2.0;
  par.T_cycle = 1.0;
  s0.P = {100.0, 95.0, 90.0, 35.0, 5.0};

  // flows temporarily zero (this will be reinitialize after reading inputs)
  s0.Q.fill(0.0);

  t_final = 20.0;
  dt = 1.0e-3;
  use_csv = true;
  inflow_csv = "AAo.csv";

  // parse file
  std::string line;
  while (std::getline(fin, line)) {
    // remove comments
    auto pos = line.find('#');
    if (pos != std::string::npos)
      line = line.substr(0, pos);

    // trim white spaces
    auto trim = [](std::string &s) {
      auto is_not_space = [](unsigned char c) { return !std::isspace(c); };
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), is_not_space));
      s.erase(std::find_if(s.rbegin(), s.rend(), is_not_space).base(), s.end());
    };
    trim(line);
    if (line.empty())
      continue;

    // expect name = value
    std::string name, value;
    std::stringstream ss(line);
    if (!std::getline(ss, name, '='))
      continue;
    if (!std::getline(ss, value))
      continue;
    trim(name);
    trim(value);

    // Numeric or string
    // ------------------ C_i ------------------
    if (!name.empty() && name[0] == 'C') {
      int idx = std::stoi(name.substr(1));
      if (0 <= idx && idx < N_NODES)
        par.C[idx] = std::stod(value);
    }
    // ------------------ R_i ------------------
    else if (!name.empty() && name[0] == 'R') {
      // R0, R1, R2, ...
      int idx = std::stoi(name.substr(1));
      if (0 <= idx && idx < N_SEGS)
        par.R[idx] = std::stod(value);
    }
    // ------------------ L_i ------------------
    else if (!name.empty() && name[0] == 'L') {
      int idx = std::stoi(name.substr(1));
      if (0 <= idx && idx < N_SEGS)
        par.L[idx] = std::stod(value);
    }
    // ------------------ P0_i ------------------
    else if (name.rfind("P0_", 0) == 0) {
      // P0_i = initial pressure at node i
      int idx = std::stoi(name.substr(3));
      if (0 <= idx && idx < N_NODES)
        s0.P[idx] = std::stod(value);
    } else if (name == "P_RA")
      par.P_RA = std::stod(value);
    else if (name == "T_cycle")
      par.T_cycle = std::stod(value);
    else if (name == "t_final")
      t_final = std::stod(value);
    else if (name == "dt")
      dt = std::stod(value);
    else if (name == "use_csv")
      use_csv = (std::stoi(value) != 0);
    else if (name == "inflow_csv")
      inflow_csv = value;
  }

  // Initialize flows after loading parameters
  for (int e = 0; e < N_SEGS; ++e) {
    int i = e;
    int j = e + 1;

    double Pup = s0.P[i];
    double Pdown = (j < N_NODES ? s0.P[j] : par.P_RA);

    if (par.L[e] > 0.0) {
      // Dynamic segment → initialize flow to 0 (safe and stable)
      s0.Q[e] = 0.0;
    } else {
      // Resistive-only segment → enforce algebraic Ohm’s law
      s0.Q[e] = (Pup - Pdown) / par.R[e];
    }
  }

  return true;
}

// 2. Inflow data aquisition
struct InflowData {
  std::vector<double> t; // time [s]
  std::vector<double> q; // flow [mL/s]
  bool loaded = false;
};

// 2-1. Default inflow function Q_in(t) [mL/s]
double Q_in_analytic(double t, const Parameters &par) {
  double omega = 2.0 * M_PI / par.T_cycle;
  return 80.0 + 40.0 * std::sin(omega * t); // mean 80 mL/s, amplitude 40 mL/s
}

// 2-2-1. Load inflow from CSV: time [s], Q [mL/s]
bool load_inflow_csv(const std::string &filename, InflowData &inflow) {
  std::ifstream fin(filename);
  if (!fin) {
    std::cerr << "Error: cannot open inflow file: " << filename << "\n";
    return false;
  }

  inflow.t.clear();
  inflow.q.clear();

  std::string line;
  bool firstLine = true;
  while (std::getline(fin, line)) {
    if (line.empty())
      continue;
    std::stringstream ss(line);
    double tt, qq;
    char sep;

    // Try "t,q" format; if it fails and it's the first line, treat it as
    // header.
    if (!(ss >> tt)) {
      if (firstLine) {
        firstLine = false;
        continue; // skip header
      } else {
        std::cerr << "Warning: could not parse line: " << line << "\n";
        continue;
      }
    }
    if (ss.peek() == ',' || ss.peek() == ';')
      ss >> sep;
    if (!(ss >> qq)) {
      std::cerr << "Warning: could not read Q on line: " << line << "\n";
      continue;
    }

    inflow.t.push_back(tt);
    inflow.q.push_back(qq);
    firstLine = false;
  }

  if (inflow.t.size() < 2) {
    std::cerr << "Error: not enough inflow points in file.\n";
    return false;
  }

  inflow.loaded = true;
  std::cout << "Loaded " << inflow.t.size() << " inflow points from "
            << filename << "\n";
  return true;
}

// 2-2-2. Linear interpolation Q_in(t) from CSV data (periodic over one cardiac
// cycle)
double Q_in_from_csv(double t, const InflowData &inflow) {
  if (!inflow.loaded) {
    throw std::runtime_error("Inflow data not loaded.");
  }
  const auto &tt = inflow.t;
  const auto &qq = inflow.q;
  double t0 = tt.front();
  double t1 = tt.back();
  double T_cycle = t1 - t0;

  // Wrap time into [t0, t1] assuming periodicity
  double t_mod = std::fmod(t - t0, T_cycle);
  if (t_mod < 0.0)
    t_mod += T_cycle;
  t_mod += t0;

  // Find interval [tt[i], tt[i+1]] containing t_mod
  // Note: Simple linear search is fine for small data; but this can be improved
  // to binary search.
  std::size_t i = 0;
  while (i + 1 < tt.size() && tt[i + 1] < t_mod) {
    ++i;
  }
  if (i + 1 >= tt.size()) {
    // Just clamp at the last point
    return qq.back();
  }

  double tL = tt[i];
  double tR = tt[i + 1];
  double qL = qq[i];
  double qR = qq[i + 1];
  double alpha = (t_mod - tL) / (tR - tL);
  return (1.0 - alpha) * qL + alpha * qR;
}

// 2.3 Selector for Q_in(t):
double Q_in(double t, const Parameters &par, const InflowData &inflow,
            bool use_csv) {
  if (use_csv && inflow.loaded) {
    return Q_in_from_csv(t, inflow);
  } else {
    return Q_in_analytic(t, par);
  }
}

// 3. Compute RHS dP/dt and dQ/dt for general node-based RLC chain
void rhs_full(double t, const State &s, const Parameters &par,
              const InflowData &inflow, bool use_csv, State &ddt) {
  double Qin = Q_in(t, par, inflow, use_csv); // [mL/s]

  // Effective flows used in mass balance (mix of dynamic & algebraic)
  std::array<double, N_SEGS> Qeff;

  // 1) dQ/dt for each segment and Qeff[e] for mass balance
  for (int e = 0; e < N_SEGS; ++e) {
    double Pup, Pdown;

    if (e < N_NODES - 1) {
      // segment e: node e -> node e+1
      Pup = s.P[e];
      Pdown = s.P[e + 1];
    } else {
      // last segment: node N_NODES-1 -> RA
      Pup = s.P[N_NODES - 1];
      Pdown = par.P_RA;
    }

    double R = par.R[e];
    double L = par.L[e];

    if (L > 0.0) {
      // Dynamic R–L segment: ODE for Q_e
      //   L_e dQ_e/dt = (P_up - P_down) - R_e Q_e
      ddt.Q[e] = ((Pup - Pdown) - R * s.Q[e]) / L;
      Qeff[e] = s.Q[e]; // for mass conservation, use dynamic Q
    } else {
      // Purely resistive segment: algebraic Ohm's law
      //   Q_e = (P_up - P_down) / R_e  (no ODE)
      ddt.Q[e] = 0.0;              // we don't integrate Q_e
      Qeff[e] = (Pup - Pdown) / R; // used in node mass balance
    }
  }

  // 2) Node mass balance using Qeff
  //    netQ[i] = sum(inflows) - sum(outflows) at node i
  std::array<double, N_NODES> netQ;
  netQ.fill(0.0);

  // Node 0: Qin enters, Q0 leaves
  netQ[0] = Qin - Qeff[0];

  // Internal nodes: i = 1..N_NODES-2
  for (int i = 1; i < N_NODES - 1; ++i) {
    netQ[i] = Qeff[i - 1] - Qeff[i];
  }

  // Last node: N_NODES-1, inflow from previous segment, outflow via last
  // segment
  netQ[N_NODES - 1] = Qeff[N_NODES - 2] - Qeff[N_SEGS - 1];

  // 3) dP/dt = netQ / C
  for (int i = 0; i < N_NODES; ++i) {
    ddt.P[i] = netQ[i] / par.C[i]; // [mmHg/s]
  }
}

int main() {
  Parameters par;
  State s, k1, k2;
  double t_final, dt;
  bool use_csv;
  std::string inflow_csv_name;

  // 1. Read input.inp for theta, ICs, and time-stepping
  if (!load_input_file("./input/input.inp", par, s, t_final, dt, use_csv,
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
    for (int e = 0; e < N_SEGS; ++e) {
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
  for (int i = 0; i < N_NODES; ++i)
    fout << "," << s.P[i];
  for (int e = 0; e < N_SEGS; ++e)
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
    for (int i = 0; i < N_NODES; ++i) {
      s_pred.P[i] = s.P[i] + dt * k1.P[i]; // predicted pressures
    }
    for (int e = 0; e < N_SEGS; ++e) {
      s_pred.Q[e] = s.Q[e] + dt * k1.Q[e]; // predicted flows
    }

    // Stage 2: evaluate RHS at predicted state (t + dt, s_pred)
    // k2 ~= dy/dt at the end of the step
    rhs_full(t + dt, s_pred, par, inflow, use_csv, k2);

    // Corrector step: average slopes k1 and k2
    // s^{n+1} = s^n + (dt/2) * (k1 + k2)
    for (int i = 0; i < N_NODES; ++i) {
      s.P[i] += 0.5 * dt * (k1.P[i] + k2.P[i]);
    }
    for (int e = 0; e < N_SEGS; ++e) {
      s.Q[e] += 0.5 * dt * (k1.Q[e] + k2.Q[e]);
    }

    // advance physical time
    t += dt;

    // Compute flows for output/diagnostics
    auto Qout = compute_Qout(s, par);

    // output.csv header: time, pressures, flows for this step
    fout << t;
    for (int i = 0; i < N_NODES; ++i)
      fout << "," << s.P[i];
    for (int e = 0; e < N_SEGS; ++e)
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