#include <chrono>

#include "main_workflow.h"

int main()
{
    auto t00 = std::chrono::high_resolution_clock::now();

    // 中文: main 仅做流程调度，具体逻辑拆分到 main_*.cpp。
    // English: main only orchestrates flow; detailed logic lives in main_*.cpp.
    mainflow::SimulationContext ctx;
    if (!mainflow::LoadAndValidateInputs(ctx))
        return 1;

    if (!mainflow::InitializeSimulationState(ctx))
        return 1;
    ctx.segment_t_start = ctx.data.t;
    auto t01 = std::chrono::high_resolution_clock::now();

    mainflow::RuntimeConfig cfg = mainflow::BuildRuntimeConfig(ctx);
    mainflow::AdvanceOutputFilesForContinue(ctx, cfg);

    if (!mainflow::InitOutput(ctx))
        return 1;

    // 中文: 先写出当前段起点（新算或续算）快照。
    // English: Emit starting snapshot of this segment (fresh or resumed).
    ctx.output.WriteData(ctx.data, ctx.pr1);
    ctx.output.WriteEnergy(ctx.data, ctx.pr1);
    ctx.output.Flush();

    if (!mainflow::RunEvolutionProgram(ctx, cfg))
        return 1;

    ctx.output.Flush();
    ctx.output.Close();

    auto t02 = std::chrono::high_resolution_clock::now();
    mainflow::PrintTiming(t00, t01, t02, ctx);
    return 0;
}
