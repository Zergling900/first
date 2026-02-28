#include "3.h"
#include "void.h"

namespace
{
enum class ForceModel
{
    BeW_ABOP_LCL = 0,
};

ForceModel g_force_model = ForceModel::BeW_ABOP_LCL;
}

bool SetForceModelByName(const std::string &name)
{
    if (name == "BeW_ABOP_LCL" || name == "bew_abop_lcl" || name == "ABOP_LCL" || name == "abop_lcl")
    {
        g_force_model = ForceModel::BeW_ABOP_LCL;
        return true;
    }
    return false;
}

const char *GetForceModelName()
{
    switch (g_force_model)
    {
    case ForceModel::BeW_ABOP_LCL:
        return "BeW_ABOP_LCL";
    default:
        return "UNKNOWN";
    }
}

void Force_BeW_ABOP_LCL(const parameter1 &pr1,
                        const parameter2 &pr2_WW,
                        const parameter2 &pr2_BB,
                        const parameter2 &pr2_WB,
                        Data &data, Cell_List &cl, vector<double> &U_atom)
{
    lcl1(data, cl, pr1, pr2_WW, U_atom);
    lcl2(data, cl, pr1, pr2_WW, pr2_BB, pr2_WB, U_atom);
}

void Force_Current(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data, Cell_List &cl, vector<double> &U_atom)
{
    // Default force model for this project version.
    // Replace this single call to swap in a new potential/LCL implementation.
    switch (g_force_model)
    {
    case ForceModel::BeW_ABOP_LCL:
    default:
        Force_BeW_ABOP_LCL(pr1, pr2_WW, pr2_BB, pr2_WB, data, cl, U_atom);
        break;
    }
}
