#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "3.h"
#include "void.h"

namespace
{
struct GroupEnergyView
{
    std::string name;
    std::vector<int> ids;
};

struct GroupObs
{
    int n = 0;
    double U = 0.0;
    double K = 0.0;
    Matrix31 F = Matrix31(0.0, 0.0, 0.0);
    Matrix31 P = Matrix31(0.0, 0.0, 0.0);
    Matrix31 MV = Matrix31(0.0, 0.0, 0.0); // sum(m*v) = sum(p/s)
    double mass_sum = 0.0;
};

static Matrix31 Zero31()
{
    return Matrix31(0.0, 0.0, 0.0);
}

static double Norm31(const Matrix31 &v)
{
    return std::sqrt(v.a00 * v.a00 + v.a10 * v.a10 + v.a20 * v.a20);
}

static double AtomMass(const Atom &a, const parameter1 &pr1)
{
    if (a.name == "W")
        return pr1.mw;
    if (a.name == "Be")
        return pr1.mb;
    return 1.0;
}

static void CollectSpeciesOrder(const Data &data, std::vector<std::string> &species_order)
{
    species_order.clear();
    for (int i = 0; i < data.n; ++i)
    {
        const std::string &name = data.atoms[i].name;
        if (name.empty())
            continue;
        if (std::find(species_order.begin(), species_order.end(), name) == species_order.end())
            species_order.push_back(name);
    }
}

static void BuildGroups(const Data &data,
                        const std::vector<std::string> &species_order,
                        std::vector<GroupEnergyView> &groups)
{
    groups.clear();

    GroupEnergyView all;
    all.name = "ALL";
    all.ids.resize(data.n);
    for (int i = 0; i < data.n; ++i)
        all.ids[i] = i;
    groups.push_back(all);

    for (size_t s = 0; s < species_order.size(); ++s)
    {
        GroupEnergyView g;
        g.name = species_order[s];
        for (int i = 0; i < data.n; ++i)
        {
            if (data.atoms[i].name == species_order[s])
                g.ids.push_back(i);
        }
        groups.push_back(g);
    }
}

static GroupObs AccumulateGroupObs(const Data &data, const parameter1 &pr1, const GroupEnergyView &g)
{
    GroupObs out;
    out.n = static_cast<int>(g.ids.size());

    const double s0 = (data.s0 == 0.0) ? 1.0 : data.s0;
    const double s02 = s0 * s0;

    for (size_t ii = 0; ii < g.ids.size(); ++ii)
    {
        const int i = g.ids[ii];
        const Atom &a = data.atoms[i];
        const double m = AtomMass(a, pr1);

        if (i >= 0 && i < static_cast<int>(data.U_atom_last.size()))
            out.U += data.U_atom_last[i];

        const double p2 = a.p.a00 * a.p.a00 + a.p.a10 * a.p.a10 + a.p.a20 * a.p.a20;
        out.K += p2 / (2.0 * s02 * m);
        out.F = out.F + a.f;
        out.P = out.P + a.p;
        out.MV = out.MV + (a.p * (1.0 / s0));
        out.mass_sum += m;
    }

    return out;
}

static void WriteEnergyHeader(FILE *fp)
{
    // time step n E H T_inst T0 s ps are kept first for easy resume parsing.
    fprintf(fp,
            "# %10s %12s %8s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n",
            "time_fs", "step", "n",
            "E_eV", "H_eV", "T_inst_K", "T0_K", "s", "ps",
            "xi", "Q", "U_eV", "K_eV",
            "F_x", "F_y", "F_z", "F_norm",
            "P_x", "P_y", "P_z", "P_norm",
            "Vcm_x", "Vcm_y", "Vcm_z");
}

static void WriteDataFrame(FILE *fp1, const Data &data, const parameter1 &pr1)
{
    double mw1 = 1.0 / pr1.mw;
    double mb1 = 1.0 / pr1.mb;
    double s0 = data.s0;
    if (s0 == 0.0)
        s0 = 1.0;

    double E = data.E;
    double E_0 = pr1.E0;

    fprintf(fp1, "%d\n", data.n);
    E = E * E_0;
    fprintf(fp1, "   time=   %8f (fs)  Energy=  %8f (eV)\n", data.t, E);
    Matrix33 Box = data.Box;

    fprintf(fp1, "BOX %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            Box.a00, Box.a10, Box.a20,
            Box.a01, Box.a11, Box.a21,
            Box.a02, Box.a12, Box.a22);

    Matrix31 v, dv;
    for (int i = 0; i < data.n; i++)
    {
        string name = data.atoms[i].name;
        Matrix31 r = data.atoms[i].r;
        Matrix31 p = data.atoms[i].p;
        Matrix31 f = data.atoms[i].f;

        if (name == "W")
        {
            v = p * mw1 * (1.0 / s0);
            dv = f * mw1;
        }
        else if (name == "Be")
        {
            v = p * mb1 * (1.0 / s0);
            dv = f * mb1;
        }
        else
        {
            v = p * (1.0 / s0);
            dv = f;
        }

        fprintf(fp1, " %4s  %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %4d\n",
                name.c_str(), r.a00, r.a10, r.a20, v.a00, v.a10, v.a20, dv.a00, dv.a10, dv.a20, data.atoms[i].identity);
    }
}

static void WriteEnergyRow(FILE *fp,
                           const GroupObs &gobs,
                           const Data &data,
                           const parameter1 &pr1,
                           bool use_global_temperature)
{
    const double E_0 = pr1.E0;
    const double s0 = data.s0;
    const double ps0 = data.ps0;
    const double t_fs = data.t;
    const long long step_id = (pr1.dt > 0.0) ? static_cast<long long>(std::llround(data.t / pr1.dt)) : 0LL;

    const double U_eV = gobs.U * E_0;
    const double K_eV = gobs.K * E_0;
    const double E_eV = (gobs.U + gobs.K) * E_0;
    const double H_eV = data.H * E_0;

    double T_inst = 0.0;
    if (use_global_temperature)
    {
        T_inst = data.T;
    }
    else if (gobs.n > 0 && pr1.kb > 0.0)
    {
        T_inst = 2.0 * gobs.K / (3.0 * gobs.n * pr1.kb);
    }

    Matrix31 Vcm = Zero31();
    if (gobs.mass_sum > 0.0)
        Vcm = gobs.MV * (1.0 / gobs.mass_sum);

    fprintf(fp,
            "%10.4f %12lld %8d %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            t_fs, step_id, gobs.n,
            E_eV, H_eV, T_inst, pr1.T, s0, ps0,
            pr1.xi, pr1.Q, U_eV, K_eV,
            gobs.F.a00, gobs.F.a10, gobs.F.a20, Norm31(gobs.F),
            gobs.P.a00, gobs.P.a10, gobs.P.a20, Norm31(gobs.P),
            Vcm.a00, Vcm.a10, Vcm.a20);
}

static bool OpenFileBuffered(const std::string &path, FILE *&fp)
{
    fp = fopen(path.c_str(), "w");
    if (!fp)
    {
        std::cerr << "Cannot open file: " << path << "\n";
        return false;
    }
    std::setvbuf(fp, nullptr, _IOFBF, 1 << 20);
    return true;
}

static bool OpenGroupedEtFiles(const std::string &base_et_path,
                               const std::vector<std::string> &species_order,
                               std::vector<FILE *> &fps,
                               std::vector<std::string> &names)
{
    fps.clear();
    names.clear();

    names.push_back("ALL");
    for (size_t i = 0; i < species_order.size(); ++i)
        names.push_back(species_order[i]);

    fps.resize(names.size(), nullptr);
    for (size_t gi = 0; gi < names.size(); ++gi)
    {
        const std::string path = MakeEtGroupPath(base_et_path, static_cast<int>(gi));
        if (!OpenFileBuffered(path, fps[gi]))
        {
            for (size_t j = 0; j < fps.size(); ++j)
            {
                if (fps[j])
                    fclose(fps[j]);
            }
            fps.clear();
            names.clear();
            return false;
        }
        fprintf(fps[gi], "# group_index=%d group_name=%s\n", static_cast<int>(gi), names[gi].c_str());
        WriteEnergyHeader(fps[gi]);
    }
    return true;
}

static void CloseGroupedEtFiles(std::vector<FILE *> &fps)
{
    for (size_t i = 0; i < fps.size(); ++i)
    {
        if (fps[i])
        {
            fclose(fps[i]);
            fps[i] = nullptr;
        }
    }
    fps.clear();
}

static void FlushGroupedEtFiles(std::vector<FILE *> &fps)
{
    for (size_t i = 0; i < fps.size(); ++i)
    {
        if (fps[i])
            fflush(fps[i]);
    }
}

static void WriteGroupedEnergyRows(std::vector<FILE *> &fps,
                                   const Data &data,
                                   const parameter1 &pr1,
                                   const std::vector<std::string> &species_order)
{
    if (fps.empty())
        return;

    std::vector<GroupEnergyView> groups;
    BuildGroups(data, species_order, groups);
    if (groups.size() != fps.size())
    {
        std::cerr << "[Output] group count changed during run; ET group files not rewritten.\n";
        return;
    }

    for (size_t gi = 0; gi < groups.size(); ++gi)
    {
        GroupObs obs = AccumulateGroupObs(data, pr1, groups[gi]);
        WriteEnergyRow(fps[gi], obs, data, pr1, gi == 0);
    }
}

} // namespace

std::string MakeEtGroupPath(const std::string &et_path, int group_index)
{
    const std::string tag = ".Et.";
    std::size_t pos = et_path.rfind(tag);
    if (pos != std::string::npos)
    {
        std::ostringstream oss;
        oss << et_path.substr(0, pos) << ".Et" << group_index << et_path.substr(pos + 3);
        return oss.str();
    }

    const std::string tag2 = ".et.";
    pos = et_path.rfind(tag2);
    if (pos != std::string::npos)
    {
        std::ostringstream oss;
        oss << et_path.substr(0, pos) << ".et" << group_index << et_path.substr(pos + 3);
        return oss.str();
    }

    std::size_t dot = et_path.find_last_of('.');
    std::ostringstream oss;
    if (dot == std::string::npos)
        oss << et_path << ".Et" << group_index;
    else
        oss << et_path.substr(0, dot) << ".Et" << group_index << et_path.substr(dot);
    return oss.str();
}

void InitEnergyFile(const FileName &filename)
{
    const std::string et0 = MakeEtGroupPath(filename.Et_file, 0);
    FILE *fp = fopen(et0.c_str(), "w");
    if (!fp)
    {
        std::cerr << "Cannot open Et0 file: " << et0 << "\n";
        return;
    }
    fprintf(fp, "# group_index=0 group_name=ALL\n");
    WriteEnergyHeader(fp);
    fclose(fp);
}

void OutputData(const Data &data, const FileName &filename, const parameter1 &pr1)
{
    FILE *fp1 = fopen(filename.Data_filename.c_str(), "a+");
    if (!fp1)
    {
        std::cerr << "Cannot open Data file: " << filename.Data_filename << "\n";
        return;
    }
    WriteDataFrame(fp1, data, pr1);
    fclose(fp1);
}

void OutputEnergy(const Data &data, const FileName &filename, const parameter1 &pr1)
{
    std::vector<std::string> species_order;
    CollectSpeciesOrder(data, species_order);

    std::vector<FILE *> fps;
    std::vector<std::string> names;
    if (!OpenGroupedEtFiles(filename.Et_file, species_order, fps, names))
        return;

    WriteGroupedEnergyRows(fps, data, pr1, species_order);
    CloseGroupedEtFiles(fps);
}

OutputManager::OutputManager()
    : data_fp_(nullptr),
      et_base_path_(),
      et_groups_initialized_(false)
{
}

OutputManager::~OutputManager()
{
    Close();
}

bool OutputManager::Open(const FileName &filename)
{
    Close();

    const std::string data_path = filename.Data_filename;
    if (!OpenFileBuffered(data_path, data_fp_))
        return false;

    et_base_path_ = filename.Et_file;
    et_group_fp_.clear();
    et_group_name_.clear();
    species_name_cache_.clear();
    et_groups_initialized_ = false;
    return true;
}

bool OutputManager::Rotate(const FileName &filename)
{
    return Open(filename);
}

void OutputManager::WriteData(const Data &data, const parameter1 &pr1)
{
    if (!data_fp_)
        return;
    WriteDataFrame(data_fp_, data, pr1);
}

void OutputManager::WriteEnergy(const Data &data, const parameter1 &pr1)
{
    if (et_base_path_.empty())
        return;

    if (!et_groups_initialized_)
    {
        CollectSpeciesOrder(data, species_name_cache_);
        if (!OpenGroupedEtFiles(et_base_path_, species_name_cache_, et_group_fp_, et_group_name_))
            return;
        et_groups_initialized_ = true;
    }

    WriteGroupedEnergyRows(et_group_fp_, data, pr1, species_name_cache_);
}

void OutputManager::Flush()
{
    if (data_fp_)
        fflush(data_fp_);
    FlushGroupedEtFiles(et_group_fp_);
}

void OutputManager::Close()
{
    if (data_fp_)
    {
        fclose(data_fp_);
        data_fp_ = nullptr;
    }

    CloseGroupedEtFiles(et_group_fp_);
    et_group_name_.clear();
    species_name_cache_.clear();
    et_base_path_.clear();
    et_groups_initialized_ = false;
}
