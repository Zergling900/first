#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "main_workflow.h"

namespace
{
bool SwitchOutputFiles(mainflow::SimulationContext &ctx, mainflow::RuntimeConfig &cfg, const std::string &reason)
{
    ++cfg.file_index;
    ctx.filename.Data_filename = mainflow::MakeIndexedPath(ctx.filename.Data_filename, cfg.file_index, cfg.data_width);
    ctx.filename.Et_file = mainflow::MakeIndexedPath(ctx.filename.Et_file, cfg.file_index, cfg.et_width);

    if (!ctx.output.Rotate(ctx.filename))
    {
        std::cerr << "Failed to rotate output files.\n";
        return false;
    }

    std::cout << "output switch (" << reason << "): Data=" << ctx.filename.Data_filename
              << " Et=" << ctx.filename.Et_file << "\n";
    return true;
}
} // namespace

namespace mainflow
{
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
    // Keep correctness unchanged while reducing I/O stall frequency.
    // If output cadence is dense, flush less often by default.
    if (cfg.output_Data_step <= 2 || cfg.output_Et_step <= 2)
        cfg.flush_every_writes = std::max(cfg.flush_every_writes, 1024);

    cfg.total_steps = EstimateTotalSteps(
        ctx.pr1,
        ctx.pr1.T,
        cfg.T_step,
        cfg.plateau_blocks,
        cfg.T_step_t,
        cfg.extra_steps_max);

    return cfg;
}

void AdvanceOutputFilesForContinue(SimulationContext &ctx, RuntimeConfig &cfg)
{
    if (!ctx.pr4.continue_run || !ctx.pr4.continue_write_next_index)
        return;

    ++cfg.file_index;
    ctx.filename.Data_filename = MakeIndexedPath(ctx.filename.Data_filename, cfg.file_index, cfg.data_width);
    ctx.filename.Et_file = MakeIndexedPath(ctx.filename.Et_file, cfg.file_index, cfg.et_width);

    std::cout << "continue mode output -> next index: Data=" << ctx.filename.Data_filename
              << " Et(base)=" << ctx.filename.Et_file
              << " Et0=" << MakeEtGroupPath(ctx.filename.Et_file, 0) << "\n";
}

// 中文: 主时间推进循环，包含注入、积分、输出与温度调度。
// English: Main time-stepping loop with injection, integration, output, and temperature schedule.
bool RunEvolutionProgram(SimulationContext &ctx, RuntimeConfig &cfg)
{
    std::cout << "STEP 10: evolution\n";

    int step_out = 0;
    int pending_writes = 0;

    auto run_steps = [&](int steps)
    {
        for (int j = 0; j < steps; ++j)
        {
            bool injected = MaybeInjectProjectiles(ctx.pr5, ctx.pr1, ctx.data, ctx.inj, false);
            if (injected)
            {
                ctx.pr1.g = std::max(0, 3 * ctx.data.n - 3);
                lcl0(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
                // 中文: 注入新粒子后先刷新受力，保证下一步半步推进一致。
                // English: Refresh forces after injection for a consistent next half-kick.
                Force_Current(ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.data, ctx.cl, ctx.U_atom);
                energy(ctx.data, ctx.pr1, ctx.U_atom);
            }

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

    const double segment_dt_fs = ctx.data.t - ctx.segment_t_start;
    if (ctx.pr1.dt > 0.0 && ms_all > 0.0)
        std::cout << "step/s = " << segment_dt_fs * 1000.0 / (ctx.pr1.dt * ms_all) << " steps/s\n";
}
} // namespace mainflow
