#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

std::string FormatIndex(int index, int width)
{
    std::ostringstream oss;
    if (width > 0)
        oss << std::setw(width) << std::setfill('0') << index;
    else
        oss << index;
    return oss.str();
}

bool AppendThermostatRange(std::ofstream &fout, int index, int width, double t_min, double t_max)
{
    if (!fout.good())
        return false;

    fout << std::fixed << std::setprecision(2)
         << t_min << "K-" << t_max << "K " << FormatIndex(index, width) << "\n";
    return fout.good();
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

    cfg.initial_hold_steps = 0;
    if (ctx.pr1.dt > 0.0 && ctx.pr3.initial_hold_time_fs > 0.0)
        cfg.initial_hold_steps = std::max(0, static_cast<int>(std::ceil(ctx.pr3.initial_hold_time_fs / ctx.pr1.dt)));

    cfg.T_step_t = std::max(1, static_cast<int>(std::ceil(ctx.pr1.T_time / ctx.pr1.dt)));
    cfg.output_Data_step = std::max(1, static_cast<int>(std::ceil(ctx.pr1.output_Data_time / ctx.pr1.dt)));
    cfg.output_Et_step = std::max(1, static_cast<int>(std::ceil(ctx.pr1.output_Et_time / ctx.pr1.dt)));
    cfg.enable_cooling = (ctx.pr3.enable_cooling != 0);
    cfg.cooling_dT = (ctx.pr3.cooling_dT > 0.0) ? ctx.pr3.cooling_dT : ctx.pr1.dT;
    cfg.cooling_T_end = (ctx.pr3.cooling_T_end >= 0.0) ? ctx.pr3.cooling_T_end : ctx.pr1.TT;
    if (cfg.cooling_dT <= 0.0)
        cfg.enable_cooling = false;

    cfg.plateau_blocks = std::max(0, ctx.pr3.plateau_blocks);
    cfg.extra_steps_max = std::max(0, ctx.pr3.extra_steps_max);
    cfg.rotate_every_data_frames = ctx.pr3.output_rotate_every_data_frames;
    cfg.flush_every_writes = std::max(1, ctx.pr3.output_flush_every_writes);
    cfg.inject_only_mode = (ctx.pr4.continue_run && ctx.pr4.resume_action == 1);
    cfg.inject_only_control_temperature = (ctx.pr5.inject_only_control_temperature != 0);
    cfg.inject_only_target_T = ctx.pr5.inject_only_target_T;
    cfg.inject_only_dT = ctx.pr5.inject_only_dT;

    double inject_only_duration_fs = ctx.pr5.inject_only_duration_fs;
    if (inject_only_duration_fs <= 0.0 && ctx.pr1.steps > 0 && ctx.pr1.dt > 0.0)
        inject_only_duration_fs = static_cast<double>(ctx.pr1.steps) * ctx.pr1.dt;
    if (cfg.inject_only_mode && inject_only_duration_fs > 0.0 && ctx.pr1.dt > 0.0)
        cfg.inject_only_steps = std::max(1, static_cast<int>(std::ceil(inject_only_duration_fs / ctx.pr1.dt)));

    // Keep correctness unchanged while reducing I/O stall frequency.
    // If output cadence is dense, flush less often by default.
    if (cfg.output_Data_step <= 2 || cfg.output_Et_step <= 2)
        cfg.flush_every_writes = std::max(cfg.flush_every_writes, 1024);

    if (cfg.inject_only_mode && cfg.inject_only_steps > 0)
    {
        cfg.total_steps = cfg.inject_only_steps;
    }
    else
    {
        cfg.total_steps = EstimateTotalSteps(
            ctx.pr1,
            ctx.pr1.T,
            cfg.T_step,
            cfg.plateau_blocks,
            cfg.T_step_t,
            cfg.extra_steps_max,
            cfg.initial_hold_steps,
            cfg.enable_cooling,
            cfg.cooling_dT,
            cfg.cooling_T_end);
    }

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

    namespace fs = std::filesystem;
    fs::path thermostat_path = fs::path(ctx.filename.Data_filename).parent_path() / "theromostat.tex";
    if (thermostat_path.empty())
        thermostat_path = fs::path("theromostat.tex");

    std::ofstream thermostat_file(thermostat_path.string(), std::ios::out | std::ios::trunc);
    if (!thermostat_file)
    {
        std::cerr << "Failed to open thermostat file: " << thermostat_path.string() << "\n";
        return false;
    }
    thermostat_file << "# T_range file_index\n";

    int step_out = 0;
    int pending_writes = 0;
    int data_frames_in_file = 1; // main.cpp writes the initial frame before this loop.
    bool range_has_data = true;
    double range_t_min = ctx.pr1.T;
    double range_t_max = ctx.pr1.T;

    auto mark_range = [&](double t_now)
    {
        if (!range_has_data)
        {
            range_t_min = t_now;
            range_t_max = t_now;
            range_has_data = true;
            return;
        }
        range_t_min = std::min(range_t_min, t_now);
        range_t_max = std::max(range_t_max, t_now);
    };

    auto flush_range = [&]() -> bool
    {
        if (!range_has_data)
            return true;
        if (!AppendThermostatRange(thermostat_file, cfg.file_index, cfg.data_width, range_t_min, range_t_max))
        {
            std::cerr << "Failed writing thermostat range for index " << cfg.file_index << "\n";
            return false;
        }
        range_has_data = false;
        return true;
    };

    auto run_steps = [&](int steps) -> bool
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

            bool wrote_data = false;
            if (step_out % cfg.output_Data_step == 0)
            {
                ctx.output.WriteData(ctx.data, ctx.pr1);
                ++pending_writes;
                ++data_frames_in_file;
                wrote_data = true;
                mark_range(ctx.pr1.T);
            }
            if (step_out % cfg.output_Et_step == 0)
            {
                ctx.output.WriteEnergy(ctx.data, ctx.pr1);
                ++pending_writes;
                mark_range(ctx.pr1.T);
            }

            if (cfg.rotate_every_data_frames > 0 && wrote_data && data_frames_in_file >= cfg.rotate_every_data_frames)
            {
                ctx.output.Flush();
                pending_writes = 0;

                if (!flush_range())
                    return false;

                std::ostringstream reason;
                reason << "every " << cfg.rotate_every_data_frames << " data frames";
                if (!SwitchOutputFiles(ctx, cfg, reason.str()))
                    return false;

                data_frames_in_file = 0;
            }

            if (pending_writes >= cfg.flush_every_writes)
            {
                ctx.output.Flush();
                pending_writes = 0;
            }
        }
        return true;
    };

    bool printed_estimate = false;

    if (cfg.inject_only_mode)
    {
        std::cout << "inject-only mode: steps=" << cfg.inject_only_steps << "\n";
        if (!cfg.inject_only_control_temperature)
        {
            ctx.pr1.xi = 0.0;
            std::cout << "inject-only mode: thermostat disabled (xi=0)\n";
        }
        else if (cfg.inject_only_target_T >= 0.0)
        {
            energy(ctx.data, ctx.pr1, ctx.U_atom);
            ctx.pr1.T = cfg.inject_only_target_T;
            ctx.pr1.H0 = ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
        }

        if (cfg.inject_only_steps <= 0)
        {
            std::cerr << "inject-only mode requested but inject_only_steps <= 0\n";
            return false;
        }

        const int temp_update_steps = std::max(1, cfg.T_step_t);
        int remaining = cfg.inject_only_steps;
        while (remaining > 0)
        {
            int block_steps = remaining;
            if (cfg.inject_only_control_temperature && cfg.inject_only_dT != 0.0)
                block_steps = std::min(block_steps, temp_update_steps);

            auto block_start = std::chrono::high_resolution_clock::now();
            if (!run_steps(block_steps))
                return false;
            auto block_end = std::chrono::high_resolution_clock::now();

            if (!printed_estimate)
            {
                const double block_ms = std::chrono::duration<double, std::milli>(block_end - block_start).count();
                double est_ms = block_ms;
                if (block_steps > 0 && cfg.total_steps > 0)
                    est_ms = block_ms * (static_cast<double>(cfg.total_steps) / static_cast<double>(block_steps));
                const double est_s = est_ms / 1000.0;
                const double est_h = est_s / 3600.0;
                std::cout << "estimate_total_time = " << est_s << " s (" << est_h << " h)\n";
                printed_estimate = true;
            }

            remaining -= block_steps;
            if (cfg.inject_only_control_temperature && cfg.inject_only_dT != 0.0)
            {
                energy(ctx.data, ctx.pr1, ctx.U_atom);
                ctx.pr1.T += cfg.inject_only_dT;
                ctx.pr1.H0 = ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
            }
        }

        if (pending_writes > 0)
            ctx.output.Flush();
        if (!flush_range())
            return false;
        thermostat_file.flush();
        return true;
    }

    if (cfg.initial_hold_steps > 0)
    {
        if (!run_steps(cfg.initial_hold_steps))
            return false;
    }

    for (int i = 0; i < cfg.T_step; ++i)
    {
        int hold_steps = HoldStepsForTemp(ctx.pr1, ctx.pr1.T, cfg.T_step_t, cfg.extra_steps_max);
        auto block_start = std::chrono::high_resolution_clock::now();

        if (!run_steps(hold_steps))
            return false;

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
        }
    }

    if (cfg.plateau_blocks > 0)
    {
        for (int p = 0; p < cfg.plateau_blocks; ++p)
        {
            int hold_steps = HoldStepsForTemp(ctx.pr1, ctx.pr1.T, cfg.T_step_t, cfg.extra_steps_max);
            if (!run_steps(hold_steps))
                return false;
        }
    }

    const double cooling_target = std::min(ctx.pr1.T, cfg.cooling_T_end);
    if (cfg.enable_cooling && cfg.cooling_dT > 0.0 && ctx.pr1.T > cooling_target + 1.0e-12)
    {
        while (ctx.pr1.T > cooling_target + 1.0e-12)
        {
            int hold_steps = HoldStepsForTemp(ctx.pr1, ctx.pr1.T, cfg.T_step_t, cfg.extra_steps_max);
            if (!run_steps(hold_steps))
                return false;

            energy(ctx.data, ctx.pr1, ctx.U_atom);
            ctx.pr1.T = std::max(ctx.pr1.T - cfg.cooling_dT, cooling_target);
            ctx.pr1.H0 = ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
        }
    }

    if (pending_writes > 0)
        ctx.output.Flush();
    if (!flush_range())
        return false;
    thermostat_file.flush();

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
