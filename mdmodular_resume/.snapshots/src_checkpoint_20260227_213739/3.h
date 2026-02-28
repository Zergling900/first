#pragma once

#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

//using namespace Eigen;
using namespace std;

// ----------------------------
//matrix
struct Matrix33
{
    double
    a00,a01,a02,
    a10,a11,a12,
    a20,a21,a22;
    Matrix33()
    : a00(0), a01(0), a02(0),
      a10(0), a11(0), a12(0),
      a20(0), a21(0), a22(0)
{}
    Matrix33(double m00, double m01, double m02,
             double m10, double m11, double m12,
             double m20, double m21, double m22)
        : a00(m00), a01(m01), a02(m02),
          a10(m10), a11(m11), a12(m12),
          a20(m20), a21(m21), a22(m22) {}
};

struct Matrix31
{
    double
    a00,a10,a20;
    Matrix31()
    : a00(0.0), a10(0.0), a20(0.0) {}

    Matrix31(double m00, double m01, double m02)
        : a00(m00), a10(m01), a20(m02) {}
};

//---------------------------------------------------------------------------------------------------------
//Matrix calculate
//33 * 33
Matrix33 operator*(const Matrix33 &A,const Matrix33 &B);
//33 +33
Matrix33 operator+(const Matrix33 &A,const Matrix33 &B);
//33 -33
Matrix33 operator-(const Matrix33 &A,const Matrix33 &B);
//33 *31
Matrix31 operator*(const Matrix33 &A,const Matrix31 &B);
//31 +31
Matrix31 operator+(const Matrix31 &A,const Matrix31 &B);
//31 - 31
Matrix31 operator-(const Matrix31 &A,const Matrix31 &B);
//k*31
Matrix31 operator*(double k, const Matrix31 &v);
//31*k
Matrix31 operator*(const Matrix31 &v, double k);
//k*33
Matrix33 operator*(double k, const Matrix33 &M);
//33*k
Matrix33 operator*(const Matrix33 &M, double k);
//---------------------------------------------------------------------------------------------------------

// struct AtomPotential
// {
    
// };
// Atom type tags for hot-path branching:
// 0 = W, 1 = Be, 2 = other/unknown.
constexpr unsigned char ATOM_TYPE_W = 0;
constexpr unsigned char ATOM_TYPE_BE = 1;
constexpr unsigned char ATOM_TYPE_OTHER = 2;

inline unsigned char AtomTypeFromName(const std::string &name)
{
    if (name == "W")
        return ATOM_TYPE_W;
    if (name == "Be")
        return ATOM_TYPE_BE;
    return ATOM_TYPE_OTHER;
}

struct Atom
{
    string name;
    int identity = 0;                // 0: original lattice, 1: injected (extensible)
    unsigned char atom_type = ATOM_TYPE_OTHER; // integer species id used in hot loops
    Matrix31 r;                     //position                    //velocity
    Matrix31 p;                     //momentum
    Matrix31 f;                     //force
};

struct Cell_List
{
    int Mx,My,Mz;
    double Wx,Wy,Wz;
    int cell_num;
    std::vector<int> Cell;
    std::vector<int> cell_offset;
    std::vector<int> atom_indices;
    // Per-cell precomputed 27-neighbor stencil (deduplicated for PBC).
    std::vector<int> cell_neighbor_offset;
    std::vector<int> cell_neighbor_cells;
    // Reusable buffers to avoid per-step allocation in lcl1.
    std::vector<int> num_in_cell_cache;
    std::vector<int> offset_cache;
    // Verlet neighbor list state.
    double verlet_skin = 0.35;
    double verlet_cutoff = 0.0;
    double verlet_cutoff2 = 0.0;
    double verlet_rebuild_limit = 0.0; // usually 0.5 * skin
    bool verlet_valid = false;
    int verlet_ref_n = 0;
    long long verlet_build_count = 0;
    std::vector<Matrix31> verlet_ref_pos;
    std::vector<int> verlet_offset;       // CSR offsets, size n+1
    std::vector<int> verlet_neighbors;    // CSR indices
    std::vector<int> verlet_count_cache;  // per-atom count scratch
    std::vector<int> verlet_cursor_cache; // per-atom cursor scratch
    // Rebuild scratch: accepted unique pairs (i<j) before CSR fill.
    std::vector<int> verlet_pair_i_cache;
    std::vector<int> verlet_pair_j_cache;
    // Persistent SoA caches for hot force kernel.
    std::vector<unsigned char> type_cache;
    std::vector<double> x_cache, y_cache, z_cache;
    std::vector<double> fx_cache, fy_cache, fz_cache, u_cache;
};

struct Data
{
    int n;                              // number of atoms
    double t,T,E,H,s0,ps0;                           // unknown
    double U_all,K_all,f_all;                           // energy
    Matrix31 F_all,P_all;
    Matrix33 Box;
    std::vector<double> U_atom_last;
    std::vector<Atom> atoms;                       //momentum
};

struct FileName
{
    string BasicData_filename,Data_filename,Et_file;
    string parameter_filename,parameter2_filename,parameter3_filename,parameter4_filename,parameter5_filename;
};

struct parameter1
{
    double dt, epsilon, kb,T, TT,T_end, sigma, mw, mb, endtime;
    double output_Data_time, output_Et_time, T_time, dT;
    double s0,ps0, xi, Q;
    double E0,H0;
    int g;
    int steps,steps_space;
};

struct parameter2
{
    double D0, r0, beta, S;
    double gamma, c, d, h;
    double R, D;
    double mu, rf, bf;
};

struct parameter3
{
    int extra_steps_max = 400000;
    int plateau_blocks = 1;
    int output_rotate_every_temp_up = 10;
    int output_rotate_every_temp_down = 10;
    int output_flush_every_writes = 1024;
    double verlet_skin = 0.50; // neighbor-list skin (distance buffer)
};

struct parameter4
{
    int continue_run = 0;              // 0: start from BasicData, 1: continue from Data/Et files
    int resume_use_target_time = 0;    // 0: use last frame/row, 1: use resume_time_fs
    double resume_time_fs = -1.0;      // target time in fs when resume_use_target_time=1
    int resume_require_exact_time = 0; // 1: fail if exact ET time not found
    int continue_write_next_index = 1; // 1: write to next xx index instead of overwriting source
    std::string force_model = "BeW_ABOP_LCL"; // runtime force dispatch selector
};

struct parameter5
{
    int enable_injection = 0;
    std::string incident_species = "W";  // exact: W / Be / He
    int incident_identity = 1;           // distinguish projectile from substrate

    double inject_start_time_fs = 0.0;
    double inject_interval_fs = 10.0;
    int inject_max_particles = 0;
    int inject_batch_size = 1;

    double inject_z_margin = 0.5;        // place particles at z = z_max - margin
    double inject_x_margin = 0.0;
    double inject_y_margin = 0.0;
    double inject_speed = 0.0;           // speed magnitude along -z

    int inject_use_fixed_seed = 1;
    unsigned int inject_seed = 123456789u;
};

struct InjectionRuntime
{
    bool initialized = false;
    int injected_total = 0;
    long long inject_events = 0;
    double next_inject_time_fs = 0.0;
    std::mt19937 rng{};
};

struct ResumeState
{
    bool valid = false;
    bool exact_time_match = false;
    double source_time_fs = 0.0;
    double target_T = 0.0; // thermostat target temperature pr1.T
    double s0 = 1.0;
    double ps0 = 0.0;
};
