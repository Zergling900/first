#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

#include "3.h"
#include "void.h"

namespace
{
double AtomMassByType(unsigned char atom_type, const parameter1 &pr1)
{
    if (atom_type == ATOM_TYPE_W)
        return pr1.mw;
    if (atom_type == ATOM_TYPE_BE)
        return pr1.mb;
    return -1.0;
}

double ClampMargin(double m, double L)
{
    if (m < 0.0)
        return 0.0;
    if (m > 0.5 * L)
        return 0.5 * L;
    return m;
}

double SpeedFromEnergyKeV(double energy_keV, double mass_amu, double E0_eV_per_internal)
{
    if (energy_keV <= 0.0 || mass_amu <= 0.0 || E0_eV_per_internal <= 0.0)
        return 0.0;

    // E(eV) = E_internal * E0, and E_internal = 0.5 * m * v^2
    const double energy_eV = energy_keV * 1000.0;
    const double energy_internal = energy_eV / E0_eV_per_internal;
    return std::sqrt(2.0 * energy_internal / mass_amu);
}
} // namespace

int ProjectileSpawnCount_Current(const parameter5 &p5, const InjectionRuntime &ir, const Data &data)
{
    (void)data;
    int remaining = p5.inject_max_particles - ir.injected_total;
    if (remaining <= 0)
        return 0;
    int batch = p5.inject_batch_size;
    if (batch < 1)
        batch = 1;
    if (batch > remaining)
        batch = remaining;
    return batch;
}

Matrix31 ProjectileVelocity_Current(const parameter5 &p5, const parameter1 &pr1, const Data &data, const Atom &a_template)
{
    (void)data;

    const double mass = AtomMassByType(a_template.atom_type, pr1);
    double v = 0.0;
    if (p5.inject_energy_keV > 0.0)
    {
        v = SpeedFromEnergyKeV(p5.inject_energy_keV, mass, pr1.E0);
    }
    else if (p5.inject_speed >= 0.0)
    {
        v = std::fabs(p5.inject_speed); // legacy fallback
    }

    return Matrix31(0.0, 0.0, -v);
}

Matrix31 ProjectilePosition_Current(const parameter5 &p5, InjectionRuntime &ir, const Data &data, int local_index)
{
    (void)local_index;
    const double Lx = data.Box.a00;
    const double Ly = data.Box.a11;
    const double Lz = data.Box.a22;

    const double mx = ClampMargin(p5.inject_x_margin, Lx);
    const double my = ClampMargin(p5.inject_y_margin, Ly);
    const double z_margin = (p5.inject_z_margin < 0.0) ? 0.0 : p5.inject_z_margin;

    double xmin = mx;
    double xmax = Lx - mx;
    if (xmax < xmin)
        xmax = xmin;

    double ymin = my;
    double ymax = Ly - my;
    if (ymax < ymin)
        ymax = ymin;

    std::uniform_real_distribution<double> ux(xmin, xmax);
    std::uniform_real_distribution<double> uy(ymin, ymax);

    double z = Lz - z_margin;
    if (z >= Lz)
        z = std::nextafter(Lz, 0.0);
    if (z < 0.0)
        z = 0.0;

    return Matrix31(ux(ir.rng), uy(ir.rng), z);
}

bool MaybeInjectProjectiles(const parameter5 &p5,
                           const parameter1 &pr1,
                           Data &data,
                           InjectionRuntime &ir,
                           bool verbose)
{
    if (!p5.enable_injection)
        return false;

    if (!(p5.incident_species == "W" || p5.incident_species == "Be" || p5.incident_species == "He"))
    {
        static bool warned_invalid_species = false;
        if (!warned_invalid_species)
        {
            warned_invalid_species = true;
            std::cerr << "[Inject] incident_species must be W / Be / He, got " << p5.incident_species << "\n";
        }
        return false;
    }

    if (!ir.initialized)
    {
        const bool random_mode = (p5.inject_seed < 0) || (p5.inject_use_fixed_seed == 0);
        unsigned int actual_seed = 0u;
        if (random_mode)
            actual_seed = std::random_device{}();
        else
            actual_seed = static_cast<unsigned int>(p5.inject_seed);

        ir.rng.seed(actual_seed);
        {
            std::ofstream fout("seed.txt", std::ios::out);
            if (fout)
            {
                fout << "inject_seed_config=" << p5.inject_seed << "\n";
                fout << "inject_seed_actual=" << actual_seed << "\n";
                fout << "seed_mode=" << (random_mode ? "random" : "fixed") << "\n";
                fout << "inject_use_fixed_seed=" << p5.inject_use_fixed_seed << "\n";
                fout << "incident_species=" << p5.incident_species << "\n";
                fout << "inject_start_time_fs=" << p5.inject_start_time_fs << "\n";
                fout << "inject_interval_fs=" << p5.inject_interval_fs << "\n";
                fout << "inject_max_particles=" << p5.inject_max_particles << "\n";
                fout << "inject_energy_keV=" << p5.inject_energy_keV << "\n";
                fout << "inject_speed_legacy=" << p5.inject_speed << "\n";
            }
        }

        if (p5.inject_energy_keV > 0.0 && p5.inject_speed >= 0.0)
            std::cout << "[Inject] both inject_energy_keV and legacy inject_speed are set; using inject_energy_keV.\n";
        else if (p5.inject_energy_keV <= 0.0 && p5.inject_speed >= 0.0)
            std::cout << "[Inject] using legacy inject_speed fallback; prefer inject_energy_keV.\n";

        ir.next_inject_time_fs = p5.inject_start_time_fs;
        ir.initialized = true;
    }

    if (p5.inject_max_particles <= 0)
        return false;

    if (p5.incident_species == "He")
    {
        static bool warned_he = false;
        if (!warned_he)
        {
            warned_he = true;
            std::cerr << "[Inject] He injection interface is reserved, but mass/potential support is not wired yet. "
                         "Use incident_species = W for validation run.\n";
        }
        return false;
    }

    const double interval = (p5.inject_interval_fs > 0.0) ? p5.inject_interval_fs : 0.0;
    const double eps = 1.0e-12;
    bool injected_any = false;

    while (ir.injected_total < p5.inject_max_particles && data.t + eps >= ir.next_inject_time_fs)
    {
        int spawn_n = ProjectileSpawnCount_Current(p5, ir, data);
        if (spawn_n <= 0)
        {
            ir.next_inject_time_fs += interval;
            ++ir.inject_events;
            if (interval <= 0.0)
                break;
            continue;
        }

        for (int k = 0; k < spawn_n; ++k)
        {
            Atom a;
            a.name = p5.incident_species;
            a.atom_type = AtomTypeFromName(a.name);
            a.identity = p5.incident_identity;
            a.r = ProjectilePosition_Current(p5, ir, data, k);
            a.f = Matrix31(0.0, 0.0, 0.0);

            Matrix31 v = ProjectileVelocity_Current(p5, pr1, data, a);
            const double m = AtomMassByType(a.atom_type, pr1);
            if (m <= 0.0)
            {
                std::cerr << "[Inject] Unsupported incident_species mass lookup for " << a.name << "\n";
                return injected_any;
            }

            const double s0 = (data.s0 == 0.0) ? 1.0 : data.s0;
            a.p = v * (m * s0);

            data.atoms.push_back(a);
            data.U_atom_last.push_back(0.0);
            ++data.n;
            ++ir.injected_total;
            injected_any = true;
        }

        if (verbose)
        {
            std::cout << "[Inject] t=" << data.t
                      << " fs event=" << ir.inject_events
                      << " added=" << spawn_n
                      << " total_injected=" << ir.injected_total
                      << " n=" << data.n << "\n";
        }

        ++ir.inject_events;
        if (interval <= 0.0)
            break;
        ir.next_inject_time_fs += interval;
    }

    return injected_any;
}
