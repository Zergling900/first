#include <iostream>
#include <chrono> // time
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#include "3.h"
#include "void.h"

namespace
{
struct SimulationContext
{
    FileName filename{};
    Data data{};
    Cell_List cl{};
    parameter1 pr1{};
    parameter2 pr2_WW{}, pr2_BB{}, pr2_WB{};
    parameter3 pr3{};
    vector<double> U_atom;
    OutputManager output;
};

struct RuntimeConfig
{
    int data_width = 2;
    int et_width = 2;
    int file_index = 0;

    int T_step = 0;
    int T_step_t = 1;
    int output_Data_step = 1;
    int output_Et_step = 1;

    int plateau_blocks = 1;
    int extra_steps_max = 10000;
    int rotate_every_up = 10;
    int rotate_every_down = 10;
    int flush_every_writes = 256;

    long long total_steps = 0;
};

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

bool ParseIndexedPath(const std::string &path, int &index, int &width)
{
    std::string base = path;
    const std::string ext = ".md3";
    if (base.size() >= ext.size() && base.compare(base.size() - ext.size(), ext.size(), ext) == 0)
    {
        base = base.substr(0, base.size() - ext.size());
    }

    std::size_t dot = base.find_last_of('.');
    if (dot == std::string::npos)
        return false;

    std::string num = base.substr(dot + 1);
    if (!IsAllDigits(num))
        return false;

    width = static_cast<int>(num.size());
    index = std::stoi(num);
    return true;
}

std::string MakeIndexedPath(const std::string &path, int index, int width)
{
    std::string base = path;
    std::string ext;
    const std::string md3 = ".md3";
    if (base.size() >= md3.size() && base.compare(base.size() - md3.size(), md3.size(), md3) == 0)
    {
        base = base.substr(0, base.size() - md3.size());
        ext = md3;
    }

    std::size_t dot = base.find_last_of('.');
    if (dot == std::string::npos)
        return path;

    std::ostringstream oss;
    oss << base.substr(0, dot + 1) << std::setw(width) << std::setfill('0') << index << ext;
    return oss.str();
}

double ComputeH0(const Data &data, const parameter1 &pr1, double target_T)
{
    return data.K_all + data.U_all
           + (data.ps0 * data.ps0) / (2.0 * pr1.Q)
           + pr1.g * pr1.kb * target_T * std::log(data.s0);
}

int HoldStepsForTemp(const parameter1 &pr1, double T, int base_steps, int extra_steps_max)
{
    double denom = pr1.T_end - pr1.TT;
    double ratio = 0.0;
    if (denom > 0.0)
        ratio = (T - pr1.TT) / denom;

    ratio = std::min(1.0, std::max(0.0, ratio));
    int extra_steps = static_cast<int>(std::lround(extra_steps_max * ratio * ratio));
    return base_steps + extra_steps;
}

long long EstimateTotalSteps(const parameter1 &pr1,
                             double T_start,
                             int heating_blocks,
                             int plateau_blocks,
                             int base_steps,
                             int extra_steps_max)
{
    long long total_steps = 0;
    double T_sim = T_start;

    for (int i = 0; i < heating_blocks; ++i)
    {
        total_steps += HoldStepsForTemp(pr1, T_sim, base_steps, extra_steps_max);
        if (T_sim < pr1.T_end)
            T_sim = std::min(T_sim + pr1.dT, pr1.T_end);
    }

    for (int p = 0; p < plateau_blocks; ++p)
    {
        total_steps += HoldStepsForTemp(pr1, T_sim, base_steps, extra_steps_max);
    }

    if (pr1.T_end > pr1.TT && pr1.dT > 0.0)
    {
        while (T_sim > pr1.TT)
        {
            total_steps += HoldStepsForTemp(pr1, T_sim, base_steps, extra_steps_max);
            T_sim = std::max(T_sim - pr1.dT, pr1.TT);
        }
    }

    return total_steps;
}

bool LoadAndValidateInputs(SimulationContext &ctx)
{
    std::cout << "STEP 1: readF\n";
    readF("FileName.txt", ctx.filename);

    if (ctx.filename.parameter_filename.empty() ||
        ctx.filename.parameter2_filename.empty() ||
        ctx.filename.BasicData_filename.empty() ||
        ctx.filename.Data_filename.empty() ||
        ctx.filename.Et_file.empty())
    {
        std::cerr << "Missing entries in FileName.txt or file not found. Check working directory.\n";
        return false;
    }

    std::cout << "STEP 2: read0\n";
    read0(ctx.filename, ctx.data);
    if (ctx.data.n <= 0 || ctx.data.atoms.size() != static_cast<size_t>(ctx.data.n))
    {
        std::cerr << "Basic data not loaded; check " << ctx.filename.BasicData_filename << ".\n";
        return false;
    }

    std::cout << "STEP 3: read1\n";
    read1(ctx.filename, ctx.pr1);
    if (ctx.pr1.dt <= 0.0)
    {
        std::cerr << "Invalid parameter1; check " << ctx.filename.parameter_filename << ".\n";
        return false;
    }

    std::cout << "STEP 4: read2\n";
    read2(ctx.filename, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
    if (ctx.pr2_WW.R <= 0.0 && ctx.pr2_BB.R <= 0.0 && ctx.pr2_WB.R <= 0.0)
    {
        std::cerr << "Invalid parameter2; check " << ctx.filename.parameter2_filename << ".\n";
        return false;
    }

    std::cout << "STEP 5: read3\n";
    read3(ctx.filename, ctx.pr3);

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

void InitializeSimulationState(SimulationContext &ctx)
{
    std::cout << "STEP 7: Dimensionalless\n";
    Dimensionalless(ctx.data, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);

    std::cout << "STEP 8: RandomV0\n";
    IntialData(ctx.data, ctx.pr1);

    std::cout << "STEP 9: First_potential\n";
    lcl0(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
    std::cout << "STEP 9: lcl0\n";
    lcl1(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.U_atom);
    std::cout << "STEP 9: lcl1\n";
    lcl2(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.U_atom);
    std::cout << "STEP 9: lcl2\n";

    energy(ctx.data, ctx.pr1, ctx.U_atom);
    ctx.pr1.H0 = ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
    energy(ctx.data, ctx.pr1, ctx.U_atom);

    ctx.output.WriteData(ctx.data, ctx.pr1);
    ctx.output.WriteEnergy(ctx.data, ctx.pr1);
    ctx.output.Flush();
}

RuntimeConfig BuildRuntimeConfig(const SimulationContext &ctx)
{
    RuntimeConfig cfg;

    int data_index = 0;
    ParseIndexedPath(ctx.filename.Data_filename, data_index, cfg.data_width);

    int et_index = 0;
    ParseIndexedPath(ctx.filename.Et_file, et_index, cfg.et_width);

    cfg.file_index = std::max(data_index, et_index);

    if (ctx.pr1.dT > 0.0 && ctx.pr1.T_end > ctx.pr1.TT)
    {
        cfg.T_step = static_cast<int>(std::ceil((ctx.pr1.T_end - ctx.pr1.TT) / ctx.pr1.dT));
    }

    cfg.T_step_t = std::max(1, static_cast<int>(std::ceil(ctx.pr1.T_time / ctx.pr1.dt)));
    cfg.output_Data_step = std::max(1, static_cast<int>(std::ceil(ctx.pr1.output_Data_time / ctx.pr1.dt)));
    cfg.output_Et_step = std::max(1, static_cast<int>(std::ceil(ctx.pr1.output_Et_time / ctx.pr1.dt)));

    cfg.plateau_blocks = std::max(0, ctx.pr3.plateau_blocks);
    cfg.extra_steps_max = std::max(0, ctx.pr3.extra_steps_max);
    cfg.rotate_every_up = ctx.pr3.output_rotate_every_temp_up;
    cfg.rotate_every_down = ctx.pr3.output_rotate_every_temp_down;
    cfg.flush_every_writes = std::max(1, ctx.pr3.output_flush_every_writes);

    cfg.total_steps = EstimateTotalSteps(
        ctx.pr1,
        ctx.pr1.T,
        cfg.T_step,
        cfg.plateau_blocks,
        cfg.T_step_t,
        cfg.extra_steps_max);

    return cfg;
}

bool SwitchOutputFiles(SimulationContext &ctx, RuntimeConfig &cfg, const std::string &reason)
{
    ++cfg.file_index;
    ctx.filename.Data_filename = MakeIndexedPath(ctx.filename.Data_filename, cfg.file_index, cfg.data_width);
    ctx.filename.Et_file = MakeIndexedPath(ctx.filename.Et_file, cfg.file_index, cfg.et_width);

    if (!ctx.output.Rotate(ctx.filename))
    {
        std::cerr << "Failed to rotate output files.\n";
        return false;
    }

    std::cout << "output switch (" << reason << "): Data=" << ctx.filename.Data_filename
              << " Et=" << ctx.filename.Et_file << "\n";
    return true;
}

bool RunEvolutionProgram(SimulationContext &ctx, RuntimeConfig &cfg)
{
    std::cout << "STEP 10: evolution\n";

    int step_out = 0;
    int pending_writes = 0;

    auto run_steps = [&](int steps)
    {
        for (int j = 0; j < steps; ++j)
        {
            ctx.data.t += ctx.pr1.dt;
            ++step_out;

            // BeW_evolution1(ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.data, ctx.cl, ctx.U_atom);
            BeW_evolution2(ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.data, ctx.cl, ctx.U_atom);

            if (step_out % cfg.output_Data_step == 0)
            {
                ctx.output.WriteData(ctx.data, ctx.pr1);
                ++pending_writes;
            }
            if (step_out % cfg.output_Et_step == 0)
            {
                ctx.output.WriteEnergy(ctx.data, ctx.pr1);
                ++pending_writes;
            }

            if (pending_writes >= cfg.flush_every_writes)
            {
                ctx.output.Flush();
                pending_writes = 0;
            }
        }
    };

    bool printed_estimate = false;
    int temp_increments = 0;

    for (int i = 0; i < cfg.T_step; ++i)
    {
        int hold_steps = HoldStepsForTemp(ctx.pr1, ctx.pr1.T, cfg.T_step_t, cfg.extra_steps_max);
        auto block_start = std::chrono::high_resolution_clock::now();

        run_steps(hold_steps);

        if (!printed_estimate)
        {
            auto block_end = std::chrono::high_resolution_clock::now();
            double block_ms = std::chrono::duration<double, std::milli>(block_end - block_start).count();
            long long block_steps = hold_steps;
            double est_ms = block_ms;
            if (block_steps > 0 && cfg.total_steps > 0)
                est_ms = block_ms * (static_cast<double>(cfg.total_steps) / static_cast<double>(block_steps));

            double est_s = est_ms / 1000.0;
            double est_h = est_s / 3600.0;
            std::cout << "estimate_total_time = " << est_s << " s (" << est_h << " h)\n";
            printed_estimate = true;
        }

        if (ctx.pr1.T < ctx.pr1.T_end)
        {
            energy(ctx.data, ctx.pr1, ctx.U_atom);
            ctx.pr1.T = std::min(ctx.pr1.T + ctx.pr1.dT, ctx.pr1.T_end);
            ctx.pr1.H0 = ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);

            ++temp_increments;
            if (cfg.rotate_every_up > 0 && temp_increments % cfg.rotate_every_up == 0)
            {
                std::ostringstream reason;
                reason << "every " << cfg.rotate_every_up << " temp up";
                if (!SwitchOutputFiles(ctx, cfg, reason.str()))
                    return false;
            }
        }
    }

    if (cfg.plateau_blocks > 0)
    {
        if (!SwitchOutputFiles(ctx, cfg, "plateau"))
            return false;

        for (int p = 0; p < cfg.plateau_blocks; ++p)
        {
            int hold_steps = HoldStepsForTemp(ctx.pr1, ctx.pr1.T, cfg.T_step_t, cfg.extra_steps_max);
            run_steps(hold_steps);
        }
    }

    if (ctx.pr1.T > ctx.pr1.TT && ctx.pr1.dT > 0.0)
    {
        if (!SwitchOutputFiles(ctx, cfg, "cooling"))
            return false;

        int temp_decrements = 0;
        while (ctx.pr1.T > ctx.pr1.TT)
        {
            int hold_steps = HoldStepsForTemp(ctx.pr1, ctx.pr1.T, cfg.T_step_t, cfg.extra_steps_max);
            run_steps(hold_steps);

            energy(ctx.data, ctx.pr1, ctx.U_atom);
            ctx.pr1.T = std::max(ctx.pr1.T - ctx.pr1.dT, ctx.pr1.TT);
            ctx.pr1.H0 = ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);

            ++temp_decrements;
            if (cfg.rotate_every_down > 0 && temp_decrements % cfg.rotate_every_down == 0 && ctx.pr1.T > ctx.pr1.TT)
            {
                std::ostringstream reason;
                reason << "every " << cfg.rotate_every_down << " temp down";
                if (!SwitchOutputFiles(ctx, cfg, reason.str()))
                    return false;
            }
        }
    }

    if (pending_writes > 0)
        ctx.output.Flush();

    return true;
}

void PrintTiming(const std::chrono::high_resolution_clock::time_point &t00,
                 const std::chrono::high_resolution_clock::time_point &t01,
                 const std::chrono::high_resolution_clock::time_point &t02,
                 const SimulationContext &ctx)
{
    double ms_all = std::chrono::duration<double, std::milli>(t02 - t00).count();
    std::cout << "time_all = " << ms_all << " ms\n";

    double ms0 = std::chrono::duration<double, std::milli>(t01 - t00).count();
    std::cout << "time_intial = " << ms0 << " ms\n";

    std::cout << "step/s = " << ctx.data.t * 1000.0 / (ctx.pr1.dt * ms_all) << " steps/s\n";
}
} // namespace

int main()
{
    auto t00 = std::chrono::high_resolution_clock::now();

    SimulationContext ctx;
    if (!LoadAndValidateInputs(ctx))
        return 1;

    if (!InitOutput(ctx))
        return 1;

    InitializeSimulationState(ctx);
    auto t01 = std::chrono::high_resolution_clock::now();

    RuntimeConfig cfg = BuildRuntimeConfig(ctx);
    if (!RunEvolutionProgram(ctx, cfg))
        return 1;

    ctx.output.Flush();
    ctx.output.Close();

    auto t02 = std::chrono::high_resolution_clock::now();
    PrintTiming(t00, t01, t02, ctx);
    return 0;
}
