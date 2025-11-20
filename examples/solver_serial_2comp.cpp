#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <stdexcept>

//---------------------------------------------
/* Example ODE system */

// Q_in(t) = inflow to arterial node from heart (LV) [mL/s]

// node a: Arterial node
// Mass conservation: Ca*dPa/dt = Q_in(t) - (Pa - Pv)/R_av

// node v: Venous node
// Mass conservation: Cv*dPv/dt = (Pa - Pv)/R_av - (Pv - P_RA)/R_vRA

// node RA: Right Atrium (boundary condition)
//---------------------------------------------

/* Solver for the 0D vascular circuit model */

// 1. Parameters for the simulation (theta vector)
struct Parameters {
    double Ca;      // arterial compliance [mL/mmHg]
    double Cv;      // venous compliance [mL/mmHg]
    double R_av;    // arterial -> venous resistance [mmHg-s/mL]
    double R_vRA;   // venous -> RA resistance [mmHg-s/mL]
    double P_RA;    // right atrium pressure [mmHg]
    double T_cyclel;// cardiac cycle period [s]
};

// 2. Inflow data aquisition
struct InflowData {
    std::vector<double> t;   // time [s]
    std::vector<double> q;   // flow [mL/s]
    bool loaded = false;
};

// 2-1. Example inflow function Q_in(t) [mL/s]
double Q_in_analytic(double t, const Parameters& par) {
    double omega = 2.0 * M_PI / par.T_cyclel;
    return 5.0 + 3.0 * std::sin(omega * t); // toy example [mL/s]
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

// 3. Compute RHS of dPa/dt and dPv/dt
void rhs(double t, double Pa, double Pv, 
         const Parameters& par, 
         const InflowData& inflow,
         bool use_csv,
         double& dPa_dt, double& dPv_dt)
{
    double Qin = Q_in(t, par, inflow, use_csv); // [mL/s]
    double Q_av = (Pa - Pv)/par.R_av;           // [mL/s] 
    double Q_vRA = (Pv - par.P_RA)/par.R_vRA;   // [mL/s]

    dPa_dt = (Qin - Q_av) / par.Ca;          // [mmHg/s]
    dPv_dt = (Q_av - Q_vRA) / par.Cv;        // [mmHg/s]
}

int main() {
    
    double Qav_sum = 0.0, QvRA_sum = 0.0;
    int count = 0;
    // 1. Set parameters (from parameters.inp file later)
    Parameters par;
    par.Ca = 1.5;        // [mL/mmHg]
    par.Cv = 20.0;        // [mL/mmHg]
    par.R_av = 1.2;        // [mmHg-s/mL]
    par.R_vRA = 0.12;        // [mmHg-s/mL]
    par.P_RA = 2.0;         // [mmHg]
    par.T_cyclel = 1;       // [s] (Not used in CSV inflow)

    // 2. Inflow data from CSV
    InflowData inflow;
    bool use_csv = true; // switch between CSV inflow and analytic example

    if (use_csv) {
        std::string inflow_file =
            "/home/jaehyn/NERS570/VascularCircuit0D/input/AAo.csv";
        if (!load_inflow_csv(inflow_file, inflow)) {
            std::cerr << "Falling back to analytic inflow.\n";
            use_csv = false;
        }
    }

    // 3. Initial conditions 
    double Pa = 100.0;       // [mmHg]
    double Pv = 5.0;         // [mmHg]

    // 4. Time stepping parameters
    double t_f = 1.0 * 20.0;       // simulate 5 cycles (if T_cycle ~ 1 s)
    double dt = 1.0e-3;      // [s]
    int n_steps = static_cast<int> (t_f / dt);

    std::ofstream fout("output.csv");
    fout << "t,Pa,Pv\n";

    double t = 0.0;
    fout << t << "," << Pa << "," << Pv << "\n";

    // 5. RK2 Time stepping loop
    for (int n=0; n < n_steps; ++n) {
        double dPa1, dPv1;
        rhs(t, Pa, Pv, par, inflow, use_csv, dPa1, dPv1);

        double Pa_pred = Pa + dt * dPa1;
        double Pv_pred = Pv + dt * dPv1;

        double dPa2, dPv2;
        rhs(t + dt, Pa_pred, Pv_pred, par, inflow, use_csv, dPa2, dPv2);
        Pa += 0.5 * dt * (dPa1 + dPa2);
        Pv += 0.5 * dt * (dPv1 + dPv2);
        
        t += dt;
        fout << t << "," << Pa << "," << Pv << "\n";
        
        double Q_av  = (Pa - Pv)/par.R_av;
        double Q_vRA = (Pv - par.P_RA)/par.R_vRA;

        Qav_sum  += Q_av;
        QvRA_sum += Q_vRA;

        count++;
    }
    fout.close();
    std::cout << "Simulation completed. Results written to output.csv\n";
    std::cout << "Mean Q_av  = " << Qav_sum  / count << " mL/s\n";
    std::cout << "Mean Q_vRA = " << QvRA_sum / count << " mL/s\n";
    return 0;
}