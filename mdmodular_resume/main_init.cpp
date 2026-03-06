#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "main_workflow.h"

namespace
{
bool IsAllDigits(const std::string &s)
{
    if (s.empty())
        return false;
    for (char c : s)
    {
        if (!std::isdigit(static_cast<unsigned char>(c)))
            return false;
    }
    return true;
}

bool PromptDoubleWithFallback(const std::string &prompt, double fallback_value, double &out_value)
{
    std::cout << prompt << " [default=" << fallback_value << "]: ";
    std::string line;
    if (!std::getline(std::cin, line))
    {
        out_value = fallback_value;
        return false;
    }
    if (line.empty())
    {
        out_value = fallback_value;
        return true;
    }
    try
    {
        out_value = std::stod(line);
        return true;
    }
    catch (...)
    {
        out_value = fallback_value;
        return false;
    }
}

std::string AllocateNextSliceDirectory(const std::string &prefix)
{
    namespace fs = std::filesystem;
    std::string safe_prefix = prefix.empty() ? "Data" : prefix;
    int max_index = 0;
    bool found = false;

    for (const auto &entry : fs::directory_iterator(fs::current_path()))
    {
        if (!entry.is_directory())
            continue;
        const std::string name = entry.path().filename().string();
        if (name.rfind(safe_prefix, 0) != 0)
            continue;
        const std::string suffix = name.substr(safe_prefix.size());
        if (!IsAllDigits(suffix))
            continue;
        const int idx = std::stoi(suffix);
        if (!found || idx > max_index)
        {
            max_index = idx;
            found = true;
        }
    }

    const int next = found ? (max_index + 1) : 1;
    return safe_prefix + std::to_string(next);
}

void RedirectOutputToSliceDirectory(mainflow::SimulationContext &ctx)
{
    namespace fs = std::filesystem;
    const std::string dir_name = AllocateNextSliceDirectory(ctx.pr4.slice_output_prefix);
    fs::create_directories(dir_name);

    const fs::path data_base = fs::path(ctx.filename.Data_filename).filename();
    const fs::path et_base = fs::path(ctx.filename.Et_file).filename();

    ctx.filename.Data_filename = (fs::path(dir_name) / data_base).string();
    ctx.filename.Et_file = (fs::path(dir_name) / et_base).string();
    ctx.pr4.continue_write_next_index = 0;

    std::cout << "slice mode output directory: " << dir_name << "\n";
    std::cout << "slice mode output paths: Data=" << ctx.filename.Data_filename
              << " Et=" << ctx.filename.Et_file << "\n";
}

struct IndexedCandidate
{
    int index = 0;
    std::string path;
};

struct DataSelection
{
    bool valid = false;
    int index = 0;
    double time_fs = 0.0;
    std::string path;
    Data frame;
};

struct EtSelection
{
    bool valid = false;
    int index = 0;
    double time_fs = 0.0;
    std::string path;
    ResumeState state;
};

bool CollectIndexedCandidates(const std::string &template_path, std::vector<IndexedCandidate> &out)
{
    namespace fs = std::filesystem;
    out.clear();

    const fs::path p(template_path);
    fs::path dir = p.parent_path();
    if (dir.empty())
        dir = fs::path(".");
    if (!fs::exists(dir) || !fs::is_directory(dir))
        return false;

    const std::string file = p.filename().string();
    const std::size_t dot_ext = file.rfind('.');
    if (dot_ext == std::string::npos)
        return false;
    const std::string ext = file.substr(dot_ext);
    const std::string stem = file.substr(0, dot_ext);

    const std::size_t dot_idx = stem.rfind('.');
    if (dot_idx == std::string::npos)
        return false;
    const std::string prefix = stem.substr(0, dot_idx + 1);
    const std::string sample_idx = stem.substr(dot_idx + 1);
    if (!IsAllDigits(sample_idx))
        return false;

    for (const auto &entry : fs::directory_iterator(dir))
    {
        if (!entry.is_regular_file())
            continue;
        const std::string name = entry.path().filename().string();
        if (name.size() <= prefix.size() + ext.size())
            continue;
        if (name.compare(0, prefix.size(), prefix) != 0)
            continue;
        if (name.compare(name.size() - ext.size(), ext.size(), ext) != 0)
            continue;

        const std::string mid = name.substr(prefix.size(), name.size() - prefix.size() - ext.size());
        if (!IsAllDigits(mid))
            continue;

        IndexedCandidate c;
        c.index = std::stoi(mid);
        c.path = (dir / name).string();
        out.push_back(c);
    }

    std::sort(out.begin(), out.end(), [](const IndexedCandidate &a, const IndexedCandidate &b)
              {
                  if (a.index != b.index)
                      return a.index < b.index;
                  return a.path < b.path;
              });
    return !out.empty();
}

void BuildCandidateListWithFallback(const std::string &template_path, std::vector<IndexedCandidate> &out)
{
    if (CollectIndexedCandidates(template_path, out))
        return;

    out.clear();
    IndexedCandidate c;
    int idx = 0;
    int width = 0;
    if (mainflow::ParseIndexedPath(template_path, idx, width))
        c.index = idx;
    c.path = template_path;
    out.push_back(c);
}

bool BetterByAbsTime(double t_candidate, int idx_candidate,
                     double t_best, int idx_best,
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
    if (t_candidate > t_best)
        return true;
    if (t_candidate < t_best)
        return false;
    return idx_candidate > idx_best;
}

bool BetterByLastTime(double t_candidate, int idx_candidate,
                      double t_best, int idx_best,
                      bool have_best)
{
    if (!have_best)
        return true;
    if (t_candidate > t_best)
        return true;
    if (t_candidate < t_best)
        return false;
    return idx_candidate > idx_best;
}

bool SelectDataFromCandidates(const std::vector<IndexedCandidate> &candidates,
                              double target_time_fs,
                              double dt_fs,
                              bool use_last,
                              bool require_exact,
                              DataSelection &out)
{
    out = DataSelection{};
    bool have_best = false;
    for (const auto &c : candidates)
    {
        Data frame;
        if (!ReadDataFrameForResume(c.path, target_time_fs, dt_fs, use_last, false, frame))
            continue;

        const bool take = use_last
                              ? BetterByLastTime(frame.t, c.index, out.time_fs, out.index, have_best)
                              : BetterByAbsTime(frame.t, c.index, out.time_fs, out.index, target_time_fs, have_best);
        if (!take)
            continue;

        out.valid = true;
        out.index = c.index;
        out.time_fs = frame.t;
        out.path = c.path;
        out.frame = frame;
        have_best = true;
    }

    if (!out.valid)
        return false;

    if (!use_last && require_exact)
    {
        const double tol = 0.01 * std::fabs(dt_fs);
        if (std::fabs(out.time_fs - target_time_fs) > tol)
            return false;
    }

    return true;
}

bool SelectEtFromCandidates(const std::vector<IndexedCandidate> &candidates,
                            double target_time_fs,
                            double dt_fs,
                            bool use_last,
                            bool require_exact,
                            const parameter4 &p4,
                            EtSelection &out)
{
    out = EtSelection{};
    bool have_best = false;
    parameter4 p4_relaxed = p4;
    p4_relaxed.resume_require_exact_time = 0;

    for (const auto &c : candidates)
    {
        ResumeState st{};
        if (!ReadEt0ResumeState(c.path, target_time_fs, dt_fs, use_last, p4_relaxed, st))
            continue;

        const bool take = use_last
                              ? BetterByLastTime(st.source_time_fs, c.index, out.time_fs, out.index, have_best)
                              : BetterByAbsTime(st.source_time_fs, c.index, out.time_fs, out.index, target_time_fs, have_best);
        if (!take)
            continue;

        out.valid = true;
        out.index = c.index;
        out.time_fs = st.source_time_fs;
        out.path = c.path;
        out.state = st;
        have_best = true;
    }

    if (!out.valid)
        return false;

    if (!use_last && require_exact)
    {
        const double tol = 0.01 * std::fabs(dt_fs);
        if (std::fabs(out.time_fs - target_time_fs) > tol)
            return false;
    }

    return true;
}

// 中文: 续算时 Data 文件中的速度/加速度列要还原为内部 p/f 变量。
// English: Convert resume-frame v/dv columns back to internal p/f representation.
void ConvertResumeFrameColumnsToPF(Data &data, const parameter1 &pr1)
{
    const double s0 = (data.s0 == 0.0) ? 1.0 : data.s0;
    for (int i = 0; i < data.n; ++i)
    {
        Atom &a = data.atoms[i];
        double m = 1.0;
        if (a.atom_type == ATOM_TYPE_W)
            m = pr1.mw;
        else if (a.atom_type == ATOM_TYPE_BE)
            m = pr1.mb;

        // Data *.md3 stores v and dv=f/m in p/f columns.
        // Data *.md3 在 p/f 列里存的是 v 与 dv=f/m。
        a.p = a.p * (m * s0);
        a.f = a.f * m;
    }
}

bool InitializeFreshSimulationState(mainflow::SimulationContext &ctx)
{
    std::cout << "STEP 7: Dimensionalless\n";
    Dimensionalless(ctx.data, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);

    std::cout << "STEP 8: RandomV0\n";
    InitialData(ctx.data, ctx.pr1);

    std::cout << "STEP 9: First_force\n";
    lcl0(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
    std::cout << "STEP 9: lcl0\n";
    Force_Current(ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.data, ctx.cl, ctx.U_atom);
    std::cout << "STEP 9: Force_Current\n";

    energy(ctx.data, ctx.pr1, ctx.U_atom);
    ctx.pr1.H0 = mainflow::ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
    energy(ctx.data, ctx.pr1, ctx.U_atom);
    return true;
}

bool InitializeResumeSimulationState(mainflow::SimulationContext &ctx)
{
    std::cout << "STEP 7: Dimensionalless (parameters only)\n";
    Dimensionalless(ctx.data, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);

    const bool slice_mode = (ctx.pr4.resume_action == 1);
    bool use_last = (ctx.pr4.resume_use_target_time == 0);
    double target_time = ctx.pr4.resume_time_fs;

    if (slice_mode)
    {
        if (ctx.pr4.prompt_slice_time)
        {
            const double fallback_time = (ctx.pr4.slice_time_fs >= 0.0) ? ctx.pr4.slice_time_fs : -1.0;
            double input_time = fallback_time;
            PromptDoubleWithFallback("Input slice time fs (<0 means last frame)", fallback_time, input_time);
            target_time = input_time;
        }
        else
        {
            target_time = ctx.pr4.slice_time_fs;
        }
        use_last = (target_time < 0.0);
    }

    const bool require_exact = (!slice_mode && ctx.pr4.resume_require_exact_time) ||
                               (slice_mode && !ctx.pr4.slice_use_nearest_time);

    const std::string data_template_path = ctx.filename.Data_filename;
    const std::string et_template_path = ctx.filename.Et_file;
    const std::string et0_template_path = MakeEtGroupPath(et_template_path, 0);

    std::vector<IndexedCandidate> data_candidates;
    std::vector<IndexedCandidate> et_candidates;
    BuildCandidateListWithFallback(data_template_path, data_candidates);
    BuildCandidateListWithFallback(et0_template_path, et_candidates);

    std::cout << "STEP 8: resume scan Data files from template " << data_template_path
              << " (candidates=" << data_candidates.size() << ")\n";
    DataSelection data_sel{};
    if (!SelectDataFromCandidates(data_candidates, target_time, ctx.pr1.dt, use_last, require_exact, data_sel))
    {
        std::cerr << "Failed to locate resume Data frame across indexed files. template=" << data_template_path
                  << " target_t=" << target_time << " use_last=" << use_last << " require_exact=" << require_exact << "\n";
        return false;
    }
    ctx.data = data_sel.frame;
    ctx.filename.Data_filename = data_sel.path;
    std::cout << "resume Data selected: file=" << data_sel.path
              << " index=" << data_sel.index
              << " t=" << data_sel.time_fs << "\n";

    ResumeState rst{};
    parameter4 et_query = ctx.pr4;
    et_query.resume_require_exact_time = require_exact ? 1 : 0;
    const double et_target_time = (slice_mode && !require_exact && !use_last) ? ctx.data.t : target_time;
    std::cout << "STEP 9: resume scan Et0 files from template " << et0_template_path
              << " (candidates=" << et_candidates.size() << ", target_t=" << et_target_time << ")\n";
    EtSelection et_sel{};
    if (!SelectEtFromCandidates(et_candidates, et_target_time, ctx.pr1.dt, use_last, require_exact, et_query, et_sel))
    {
        std::cerr << "Failed to locate resume Et0 row across indexed files. template=" << et0_template_path
                  << " target_t=" << et_target_time << " use_last=" << use_last << " require_exact=" << require_exact << "\n";
        return false;
    }
    rst = et_sel.state;
    std::cout << "resume Et0 selected: file=" << et_sel.path
              << " index=" << et_sel.index
              << " t=" << et_sel.time_fs << "\n";

    int et_idx = 0;
    int et_width = 0;
    if (mainflow::ParseIndexedPath(et_template_path, et_idx, et_width))
        ctx.filename.Et_file = mainflow::MakeIndexedPath(et_template_path, et_sel.index, et_width);
    else
        ctx.filename.Et_file = et_template_path;

    if (data_sel.index != et_sel.index)
    {
        std::cout << "warning: resume Data/Et selected different file indices: Data=" << data_sel.index
                  << " Et=" << et_sel.index << "\n";
    }

    if (rst.valid)
    {
        const double tol = 0.01 * std::fabs(ctx.pr1.dt);
        if (std::fabs(ctx.data.t - rst.source_time_fs) > tol)
        {
            if (require_exact)
            {
                std::cerr << "resume mismatch between Data and Et0 time: data.t=" << ctx.data.t
                          << " et0.t=" << rst.source_time_fs << " tol=" << tol << "\n";
                return false;
            }
            std::cout << "warning: nearest Data/Et0 times differ: data.t=" << ctx.data.t
                      << " et0.t=" << rst.source_time_fs << "\n";
        }
        ctx.data.s0 = rst.s0;
        ctx.data.ps0 = rst.ps0;
        ctx.pr1.T = rst.target_T;
        std::cout << "resume state: data.t=" << ctx.data.t
                  << " et0.t=" << rst.source_time_fs
                  << " targetT=" << ctx.pr1.T
                  << " s0=" << ctx.data.s0
                  << " ps0=" << ctx.data.ps0 << "\n";
    }
    else
    {
        ctx.data.s0 = ctx.pr1.s0;
        ctx.data.ps0 = ctx.pr1.ps0;
        ctx.pr1.T = ctx.pr1.TT;
        std::cout << "warning: Et0 resume state missing, fallback to parameter values.\n";
    }

    if (slice_mode)
    {
        const double source_t = ctx.data.t;
        ctx.data.t = 0.0;
        std::cout << "slice mode: reset simulation time from " << source_t << " fs to 0 fs\n";
    }

    ctx.pr1.g = 3 * ctx.data.n - 3;
    ConvertResumeFrameColumnsToPF(ctx.data, ctx.pr1);

    std::cout << "STEP 10: Rebuild Cell/Force after resume\n";
    lcl0(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
    Force_Current(ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.data, ctx.cl, ctx.U_atom);
    energy(ctx.data, ctx.pr1, ctx.U_atom);
    ctx.pr1.H0 = mainflow::ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
    energy(ctx.data, ctx.pr1, ctx.U_atom);

    if (slice_mode)
        RedirectOutputToSliceDirectory(ctx);

    return true;
}
} // namespace

namespace mainflow
{
// 中文: 读取并校验所有输入参数与数据文件。
// English: Load and validate all input parameters and data files.
bool LoadAndValidateInputs(SimulationContext &ctx)
{
    std::cout << "STEP 1: readF\n";
    readF("FileName.txt", ctx.filename);

    if (ctx.filename.parameter_filename.empty() ||
        ctx.filename.parameter2_filename.empty() ||
        ctx.filename.Data_filename.empty() ||
        ctx.filename.Et_file.empty())
    {
        std::cerr << "Missing entries in FileName.txt or file not found. Check working directory.\n";
        return false;
    }

    std::cout << "STEP 2: read1\n";
    read1(ctx.filename, ctx.pr1);
    if (ctx.pr1.dt <= 0.0)
    {
        std::cerr << "Invalid parameter1; check " << ctx.filename.parameter_filename << ".\n";
        return false;
    }

    std::cout << "STEP 3: read2\n";
    read2(ctx.filename, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
    if (ctx.pr2_WW.R <= 0.0 && ctx.pr2_BB.R <= 0.0 && ctx.pr2_WB.R <= 0.0)
    {
        std::cerr << "Invalid parameter2; check " << ctx.filename.parameter2_filename << ".\n";
        return false;
    }

    std::cout << "STEP 4: read3\n";
    read3(ctx.filename, ctx.pr3);
    if (ctx.pr3.verlet_skin > 0.0)
        ctx.cl.verlet_skin = ctx.pr3.verlet_skin;
    std::cout << "STEP 5: read4\n";
    read4(ctx.filename, ctx.pr4);
    std::cout << "STEP 5b: read5\n";
    read5(ctx.filename, ctx.pr5);
    if (!SetForceModelByName(ctx.pr4.force_model))
    {
        std::cerr << "Unknown force_model in parameter.p4: " << ctx.pr4.force_model << "\n";
        return false;
    }
    std::cout << "force_model = " << GetForceModelName() << "\n";
    if (!(ctx.pr5.incident_species == "W" || ctx.pr5.incident_species == "Be" || ctx.pr5.incident_species == "He"))
    {
        std::cerr << "incident_species must be exact W / Be / He, got: " << ctx.pr5.incident_species << "\n";
        return false;
    }

    if (!ctx.pr4.continue_run)
    {
        std::cout << "STEP 6: read0 (BasicData)\n";
        read0(ctx.filename, ctx.data);
        if (ctx.data.n <= 0 || ctx.data.atoms.size() != static_cast<size_t>(ctx.data.n))
        {
            std::cerr << "Basic data not loaded; check " << ctx.filename.BasicData_filename << ".\n";
            return false;
        }
    }
    else
    {
        std::cout << "STEP 6: continue mode enabled (data will be loaded from Data/Et files)\n";
    }

    return true;
}

bool InitOutput(SimulationContext &ctx)
{
    std::cout << "STEP 6: InitOutput\n";
    if (!ctx.output.Open(ctx.filename))
    {
        std::cerr << "Failed to initialize output files.\n";
        return false;
    }
    std::cout << "output: Data=" << ctx.filename.Data_filename
              << " Et=" << ctx.filename.Et_file << "\n";
    return true;
}

bool InitializeSimulationState(SimulationContext &ctx)
{
    if (ctx.pr4.continue_run)
        return InitializeResumeSimulationState(ctx);
    return InitializeFreshSimulationState(ctx);
}
} // namespace mainflow
