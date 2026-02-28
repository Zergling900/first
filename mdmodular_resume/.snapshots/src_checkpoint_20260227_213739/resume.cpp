#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "3.h"
#include "void.h"

namespace
{
bool ReadOneDataFrame(std::ifstream &fin, Data &frame)
{
    int n = 0;
    if (!(fin >> n))
        return false;

    std::string dummy;
    std::getline(fin, dummy);

    std::string tmp;
    double t_fs = 0.0;
    double E_eV = 0.0;
    if (!(fin >> tmp >> t_fs >> tmp >> tmp >> E_eV))
        return false;
    std::getline(fin, dummy);

    if (!(fin >> tmp))
        return false;
    if (!(fin >> frame.Box.a00 >> frame.Box.a01 >> frame.Box.a02
              >> frame.Box.a10 >> frame.Box.a11 >> frame.Box.a12
              >> frame.Box.a20 >> frame.Box.a21 >> frame.Box.a22))
        return false;
    std::getline(fin, dummy);

    frame.n = n;
    frame.t = t_fs;
    frame.E = E_eV; // eV in file; caller may recompute, so this is informational
    frame.atoms.resize(n);

    for (int i = 0; i < n; ++i)
    {
        std::string atom_line;
        if (!std::getline(fin, atom_line))
            return false;
        if (atom_line.empty())
        {
            --i;
            continue;
        }

        Atom &a = frame.atoms[i];
        a.identity = 0;
        std::stringstream ss(atom_line);
        if (!(ss >> a.name
                  >> a.r.a00 >> a.r.a10 >> a.r.a20
                  >> a.p.a00 >> a.p.a10 >> a.p.a20
                  >> a.f.a00 >> a.f.a10 >> a.f.a20))
            return false;
        if (!(ss >> a.identity))
            a.identity = 0;
        a.atom_type = AtomTypeFromName(a.name);
    }

    return true;
}

bool IsBetterTimeCandidate(double t_candidate,
                           double t_best,
                           double target_time_fs,
                           bool have_best)
{
    if (!have_best)
        return true;

    const bool cand_le = (t_candidate <= target_time_fs);
    const bool best_le = (t_best <= target_time_fs);

    if (cand_le != best_le)
        return cand_le; // prefer <= target

    if (cand_le)
        return t_candidate > t_best; // latest not exceeding target

    return t_candidate < t_best; // otherwise earliest exceeding target
}

bool IsBetterAbsTimeCandidate(double t_candidate,
                              double t_best,
                              double target_time_fs,
                              bool have_best)
{
    if (!have_best)
        return true;
    const double dc = std::fabs(t_candidate - target_time_fs);
    const double db = std::fabs(t_best - target_time_fs);
    if (dc < db)
        return true;
    if (dc > db)
        return false;
    // Tie-breaker: prefer later time.
    return t_candidate > t_best;
}

bool ReadEtRow(std::ifstream &fin,
               double &time_fs,
               long long &step_id,
               int &n,
               double &E,
               double &H,
               double &T_inst,
               double &T_target,
               double &s0,
               double &ps0)
{
    std::string line;
    while (std::getline(fin, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '#')
            continue;
        if (line.find("time") != std::string::npos && line.find("step") != std::string::npos)
            continue; // header

        std::stringstream ss(line);

        // Column order must match output.cpp (ET grouped writer).
        // time step n E H T_inst T_target s ps ...
        if (!(ss >> time_fs >> step_id >> n >> E >> H >> T_inst >> T_target >> s0 >> ps0))
            continue;

        return true;
    }
    return false;
}
} // namespace

bool ReadDataFrameForResume(const std::string &path,
                            double target_time_fs,
                            double dt_fs,
                            bool use_last_frame,
                            Data &data)
{
    std::ifstream fin(path);
    if (!fin)
    {
        std::cerr << "can't open data file for resume: " << path << "\n";
        return false;
    }

    Data frame;
    Data best;
    bool have_best = false;
    double best_t = 0.0;

    while (ReadOneDataFrame(fin, frame))
    {
        if (use_last_frame)
        {
            best = frame;
            best_t = frame.t;
            have_best = true;
            continue;
        }

        if (IsBetterAbsTimeCandidate(frame.t, best_t, target_time_fs, have_best))
        {
            best = frame;
            best_t = frame.t;
            have_best = true;
        }
    }

    if (!have_best)
    {
        std::cerr << "no frame found in data file: " << path << "\n";
        return false;
    }

    if (!use_last_frame)
    {
        const double tol = 0.01 * std::fabs(dt_fs);
        if (std::fabs(best_t - target_time_fs) > tol)
        {
            std::cerr << "resume data frame time mismatch: target=" << target_time_fs
                      << " chosen=" << best_t << " tol=" << tol << "\n";
            return false;
        }
    }

    data = best;
    return true;
}

bool ReadEt0ResumeState(const std::string &et0_path,
                        double target_time_fs,
                        double dt_fs,
                        bool use_last_row,
                        const parameter4 &p4,
                        ResumeState &st)
{
    st = ResumeState{};

    std::ifstream fin(et0_path);
    if (!fin)
    {
        std::cerr << "can't open Et0 file for resume: " << et0_path << "\n";
        return false;
    }

    bool have_best = false;
    double best_t = 0.0;
    long long best_step = 0;
    int best_n = 0;
    double best_E = 0.0, best_H = 0.0, best_T_inst = 0.0;
    double best_T_target = 0.0, best_s0 = 1.0, best_ps0 = 0.0;

    double t = 0.0;
    long long step_id = 0;
    int n = 0;
    double E = 0.0, H = 0.0, T_inst = 0.0, T_target = 0.0, s0 = 1.0, ps0 = 0.0;

    while (ReadEtRow(fin, t, step_id, n, E, H, T_inst, T_target, s0, ps0))
    {
        if (use_last_row)
        {
            have_best = true;
            best_t = t;
            best_step = step_id;
            best_n = n;
            best_E = E;
            best_H = H;
            best_T_inst = T_inst;
            best_T_target = T_target;
            best_s0 = s0;
            best_ps0 = ps0;
            continue;
        }

        if (IsBetterAbsTimeCandidate(t, best_t, target_time_fs, have_best))
        {
            have_best = true;
            best_t = t;
            best_step = step_id;
            best_n = n;
            best_E = E;
            best_H = H;
            best_T_inst = T_inst;
            best_T_target = T_target;
            best_s0 = s0;
            best_ps0 = ps0;
        }
    }

    if (!have_best)
    {
        std::cerr << "no ET row found in: " << et0_path << "\n";
        return false;
    }

    if (!use_last_row)
    {
        const double tol = 0.01 * std::fabs(dt_fs);
        if (std::fabs(best_t - target_time_fs) > tol)
        {
            std::cerr << "resume Et0 time mismatch: target=" << target_time_fs
                      << " chosen=" << best_t << " tol=" << tol << "\n";
            return false;
        }
        if (p4.resume_require_exact_time)
            st.exact_time_match = true;
    }

    st.valid = true;
    st.exact_time_match = use_last_row ? true : (std::fabs(best_t - target_time_fs) <= 0.01 * std::fabs(dt_fs));
    st.source_time_fs = best_t;
    st.target_T = best_T_target;
    st.s0 = best_s0;
    st.ps0 = best_ps0;
    (void)best_step;
    (void)best_n;
    (void)best_E;
    (void)best_H;
    (void)best_T_inst;
    return true;
}
