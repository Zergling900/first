#include <cmath>
#include <iostream>
#include <string>

#include "main_workflow.h"

namespace
{
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
    IntialData(ctx.data, ctx.pr1);

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

    const bool use_last = (ctx.pr4.resume_use_target_time == 0);
    const double target_time = ctx.pr4.resume_time_fs;

    const std::string data_path = ctx.filename.Data_filename;
    const std::string et0_path = MakeEtGroupPath(ctx.filename.Et_file, 0);

    std::cout << "STEP 8: resume read Data frame from " << data_path << "\n";
    if (!ReadDataFrameForResume(data_path, target_time, ctx.pr1.dt, use_last, ctx.data))
        return false;

    ResumeState rst{};
    std::cout << "STEP 9: resume read Et0 row from " << et0_path << "\n";
    if (!ReadEt0ResumeState(et0_path, target_time, ctx.pr1.dt, use_last, ctx.pr4, rst))
        return false;

    if (rst.valid)
    {
        const double tol = 0.01 * std::fabs(ctx.pr1.dt);
        if (std::fabs(ctx.data.t - rst.source_time_fs) > tol)
        {
            std::cerr << "resume mismatch between Data and Et0 time: data.t=" << ctx.data.t
                      << " et0.t=" << rst.source_time_fs << " tol=" << tol << "\n";
            return false;
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

    ctx.pr1.g = 3 * ctx.data.n - 3;
    ConvertResumeFrameColumnsToPF(ctx.data, ctx.pr1);

    std::cout << "STEP 10: Rebuild Cell/Force after resume\n";
    lcl0(ctx.data, ctx.cl, ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB);
    Force_Current(ctx.pr1, ctx.pr2_WW, ctx.pr2_BB, ctx.pr2_WB, ctx.data, ctx.cl, ctx.U_atom);
    energy(ctx.data, ctx.pr1, ctx.U_atom);
    ctx.pr1.H0 = mainflow::ComputeH0(ctx.data, ctx.pr1, ctx.pr1.T);
    energy(ctx.data, ctx.pr1, ctx.U_atom);
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
