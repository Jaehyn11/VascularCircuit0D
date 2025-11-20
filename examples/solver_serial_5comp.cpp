#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cstdlib>

//---------------------------------------------
// 0D Vascular circuit model: Ao -> ART -> ARTR -> CAP -> VEN -> RA
// Units:
//   Q:     mL/s
//   P:     mmHg
//   C:     mL/mmHg
//   R:     mmHg·s/mL
//   t, dt: s
//---------------------------------------------

/* Solver for the 0D vascular circuit model */

// 1-1. Parameters for the simulation (theta vector)
struct Parameters {
    // Compliances [mL/mmHg]
    double C_ao;
    double C_art;
    double C_artr;
    double C_cap;
    double C_ven;

    // Resistances [mmHg·s/mL]
    double R_ao_art;
    double R_art_artr;
    double R_artr_cap;
    double R_cap_ven;
    double R_ven_ra;

    // Right atrial pressure [mmHg]
    double P_RA;

    // Cardiac period (only if you still use analytic inflow)
    double T_cycle;
};

// 1-2. State vector - pressures (y vector)
struct State {
    double P_ao;
    double P_art;
    double P_artr;
    double P_cap;
    double P_ven;
};


// 2. Inflow data aquisition
struct InflowData {
    std::vector<double> t;   // time [s]
    std::vector<double> q;   // flow [mL/s]
    bool loaded = false;
};

// 2-1. Example inflow function Q_in(t) [mL/s]
double Q_in_analytic(double t, const Parameters& par) {
    double omega = 2.0 * M_PI / par.T_cycle;
    // toy example: mean 80 mL/s, amplitude 40 mL/s
    return 80.0 + 40.0 * std::sin(omega * t);
}

// 2-2-1. Load inflow from CSV: time [s], Q [mL/s]
bool load_inflow_csv(const std::string& filename, InflowData& inflow)
{
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
        if (line.empty()) continue;

        std::stringstream ss(line);
        double tt, qq;
        char sep;

        // Try "t,q" format; if it fails and it's the first line, treat it as header.
        if (!(ss >> tt)) {
            if (firstLine) {
                firstLine = false;
                continue; // skip header
            } else {
                std::cerr << "Warning: could not parse line: " << line << "\n";
                continue;
            }
        }
        if (ss.peek() == ',' || ss.peek() == ';') ss >> sep;
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
    std::cout << "Loaded " << inflow.t.size() << " inflow points from " << filename << "\n";
    return true;
}

// 2-2-2. Linear interpolation Q_in(t) from CSV data (periodic over one cardiac cycle)
double Q_in_from_csv(double t, const InflowData& inflow)
{
    if (!inflow.loaded) {
        throw std::runtime_error("Inflow data not loaded.");
    }

    const auto& tt = inflow.t;
    const auto& qq = inflow.q;

    double t0 = tt.front();
    double t1 = tt.back();
    double T_cycle = t1 - t0;

    // Wrap time into [t0, t1] assuming periodicity
    double t_mod = std::fmod(t - t0, T_cycle);
    if (t_mod < 0.0) t_mod += T_cycle;
    t_mod += t0;

    // Find interval [tt[i], tt[i+1]] containing t_mod
    // (simple linear search is fine for small data; can be improved to binary search)
    std::size_t i = 0;
    while (i + 1 < tt.size() && tt[i+1] < t_mod) {
        ++i;
    }
    if (i + 1 >= tt.size()) {
        // Just clamp at the last point
        return qq.back();
    }

    double tL = tt[i];
    double tR = tt[i+1];
    double qL = qq[i];
    double qR = qq[i+1];

    double alpha = (t_mod - tL) / (tR - tL);
    return (1.0 - alpha) * qL + alpha * qR;
}

// 2.3 Selector for Q_in(t):
//    if use_csv == true, use data-driven inflow, otherwise use analytic example.
double Q_in(double t, const Parameters& par,
            const InflowData& inflow, bool use_csv)
{
    if (use_csv && inflow.loaded) {
        return Q_in_from_csv(t, inflow);
    } else {
        return Q_in_analytic(t, par);
    }
}

// 3. Compute RHS dP/dt for full 5-compartment system
void rhs_full(double t,
              const State& s,
              const Parameters& par,
              const InflowData& inflow,
              bool use_csv,
              State& ddt)
{
    // Inflow from LV [mL/s]
    double Qin = Q_in(t, par, inflow, use_csv);

    // Flows between compartments [mL/s]
    double Q_ao_art   = (s.P_ao   - s.P_art)  / par.R_ao_art;
    double Q_art_artr = (s.P_art  - s.P_artr) / par.R_art_artr;
    double Q_artr_cap = (s.P_artr - s.P_cap)  / par.R_artr_cap;
    double Q_cap_ven  = (s.P_cap  - s.P_ven)  / par.R_cap_ven;
    double Q_ven_RA   = (s.P_ven  - par.P_RA) / par.R_ven_ra;

    // dP/dt [mmHg/s] = (net flow [mL/s]) / C [mL/mmHg]
    ddt.P_ao   = (Qin        - Q_ao_art)   / par.C_ao;
    ddt.P_art  = (Q_ao_art   - Q_art_artr)/ par.C_art;
    ddt.P_artr = (Q_art_artr - Q_artr_cap)/ par.C_artr;
    ddt.P_cap  = (Q_artr_cap - Q_cap_ven) / par.C_cap;
    ddt.P_ven  = (Q_cap_ven  - Q_ven_RA)  / par.C_ven;
}


int main() {
    
    // ---- 1. Set parameters (physiologic-ish guesses) ----
    Parameters par;
    par.C_ao      = 2.0;   // mL/mmHg
    par.C_art     = 3.0;
    par.C_artr    = 1.0;
    par.C_cap     = 0.5;
    par.C_ven     = 20.0;

    par.R_ao_art  = 0.5;   // mmHg·s/mL
    par.R_art_artr= 0.5;
    par.R_artr_cap= 1.0;
    par.R_cap_ven = 1.0;
    par.R_ven_ra  = 0.2;

    par.P_RA      = 2.0;   // mmHg
    par.T_cycle   = 1.0;   // s (for analytic inflow)

    // 2. Inflow data from CSV
    InflowData inflow;
    bool use_csv = true; // switch between CSV inflow and analytic example

    if (use_csv) {
        // Build path from $HOME/VascularCircuit/input/AAo.csv
        const char* home_env = std::getenv("HOME");
        std::string inflow_file;

        if (home_env) {
            std::string home(home_env);
            inflow_file = home + "/VascularCircuit/input/AAo.csv";
        } else {
            // Fallback: look in current directory
            inflow_file = "AAo.csv";
        }

        std::cout << "Trying inflow file: " << inflow_file << "\n";

        if (!load_inflow_csv(inflow_file, inflow)) {
            std::cerr << "Falling back to analytic inflow.\n";
            use_csv = false;
        }
    }

    // 3. Initial conditions
    State s, k1, k2;
    s.P_ao   = 100.0;
    s.P_art  = 95.0;
    s.P_artr = 90.0;
    s.P_cap  = 35.0;
    s.P_ven  = 5.0;

    // 4. Time stepping parameters
    double t_f   = 20.0;       // simulate 20 s (~20 cycles)
    double dt    = 1.0e-3;     // s
    int    n_steps = static_cast<int>(t_f / dt);

    std::ofstream fout("output.csv");
    fout << "t,P_ao,P_art,P_artr,P_cap,P_ven\n";

    double t = 0.0;
    fout << t << ","
         << s.P_ao   << ","
         << s.P_art  << ","
         << s.P_artr << ","
         << s.P_cap  << ","
         << s.P_ven  << "\n";

    // For mean flow diagnostics
    double Qav_sum  = 0.0;
    double QvRA_sum = 0.0;
    int    count    = 0;

    // 5. RK2 Time stepping loop
    for (int n = 0; n < n_steps; ++n) {
        rhs_full(t, s, par, inflow, use_csv, k1);

        State s_pred;
        s_pred.P_ao   = s.P_ao   + dt * k1.P_ao;
        s_pred.P_art  = s.P_art  + dt * k1.P_art;
        s_pred.P_artr = s.P_artr + dt * k1.P_artr;
        s_pred.P_cap  = s.P_cap  + dt * k1.P_cap;
        s_pred.P_ven  = s.P_ven  + dt * k1.P_ven;

        rhs_full(t + dt, s_pred, par, inflow, use_csv, k2);

        s.P_ao   += 0.5 * dt * (k1.P_ao   + k2.P_ao);
        s.P_art  += 0.5 * dt * (k1.P_art  + k2.P_art);
        s.P_artr += 0.5 * dt * (k1.P_artr + k2.P_artr);
        s.P_cap  += 0.5 * dt * (k1.P_cap  + k2.P_cap);
        s.P_ven  += 0.5 * dt * (k1.P_ven  + k2.P_ven);

        t += dt;

        fout << t << ","
             << s.P_ao   << ","
             << s.P_art  << ","
             << s.P_artr << ","
             << s.P_cap  << ","
             << s.P_ven  << "\n";

        // Flow diagnostics at current state
        double Q_ao_art   = (s.P_ao   - s.P_art)  / par.R_ao_art;
        double Q_ven_RA   = (s.P_ven  - par.P_RA) / par.R_ven_ra;
        Qav_sum  += Q_ao_art;
        QvRA_sum += Q_ven_RA;
        count++;
    }

    fout.close();
    std::cout << "Simulation completed. Results written to output.csv\n";
    std::cout << "Mean Q_ao_art = " << Qav_sum  / count << " mL/s\n";
    std::cout << "Mean Q_ven_RA = " << QvRA_sum / count << " mL/s\n";

    return 0;
}