#pragma once

#include <chrono>
#include <string>
#include <vector>

#include "void.h"

namespace mainflow
{
// 中文: 运行上下文，集中保存一次仿真所需的全局状态。
// English: Simulation context holding all state needed by one run.
struct SimulationContext
{
    FileName filename{};
    Data data{};
    Cell_List cl{};
    parameter1 pr1{};
    parameter2 pr2_WW{}, pr2_BB{}, pr2_WB{};
    parameter3 pr3{};
    parameter4 pr4{};
    parameter5 pr5{};
    InjectionRuntime inj{};
    double segment_t_start = 0.0;
    std::vector<double> U_atom;
    OutputManager output;
};

// 中文: 运行时调度参数（输出步长、温度分段、文件轮转等）。
// English: Runtime scheduling config (output cadence, temperature blocks, file rotation).
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

bool ParseIndexedPath(const std::string &path, int &index, int &width);
std::string MakeIndexedPath(const std::string &path, int index, int width);

double ComputeH0(const Data &data, const parameter1 &pr1, double target_T);
int HoldStepsForTemp(const parameter1 &pr1, double T, int base_steps, int extra_steps_max);
long long EstimateTotalSteps(const parameter1 &pr1,
                             double T_start,
                             int heating_blocks,
                             int plateau_blocks,
                             int base_steps,
                             int extra_steps_max);

bool LoadAndValidateInputs(SimulationContext &ctx);
bool InitOutput(SimulationContext &ctx);
bool InitializeSimulationState(SimulationContext &ctx);

RuntimeConfig BuildRuntimeConfig(const SimulationContext &ctx);
void AdvanceOutputFilesForContinue(SimulationContext &ctx, RuntimeConfig &cfg);
bool RunEvolutionProgram(SimulationContext &ctx, RuntimeConfig &cfg);

void PrintTiming(const std::chrono::high_resolution_clock::time_point &t00,
                 const std::chrono::high_resolution_clock::time_point &t01,
                 const std::chrono::high_resolution_clock::time_point &t02,
                 const SimulationContext &ctx);
} // namespace mainflow
