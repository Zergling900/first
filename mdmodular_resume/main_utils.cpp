#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>

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
} // namespace

namespace mainflow
{
// 中文: 解析带数字尾缀的路径（例如 BCC.09.md3）。
// English: Parse indexed path suffix (e.g., BCC.09.md3).
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

// 中文: 用给定编号重建输出路径，保留原始位宽。
// English: Build indexed output path with fixed width.
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

// 中文: 计算 Nose/NVT 使用的守恒量 H0。
// English: Compute conserved quantity H0 for Nose/NVT workflow.
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
} // namespace mainflow
