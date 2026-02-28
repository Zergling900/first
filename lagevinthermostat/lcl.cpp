#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <omp.h>

#include "3.h"
#include "void.h"

void lcl0(Data &data, Cell_List &cl, const parameter1 &pr1,
          const parameter2 &pr2_WW,
          const parameter2 &pr2_BB,
          const parameter2 &pr2_WB)
{
    int n = data.n;
    (void)pr1;

    const double cutoff_WW = pr2_WW.R + pr2_WW.D;
    const double cutoff_BB = pr2_BB.R + pr2_BB.D;
    const double cutoff_WB = pr2_WB.R + pr2_WB.D;
    double cutoff = std::max(cutoff_WW, std::max(cutoff_BB, cutoff_WB));
    if (cutoff <= 0.0 || data.Box.a00 <= 0.0 || data.Box.a11 <= 0.0 || data.Box.a22 <= 0.0)
    {
        std::cerr << "[lcl0] Invalid box or cutoff (R+D).\n";
        cl = Cell_List{};
        return;
    }

    cl.Mx = static_cast<int>(std::floor(data.Box.a00 / cutoff));
    cl.My = static_cast<int>(std::floor(data.Box.a11 / cutoff));
    cl.Mz = static_cast<int>(std::floor(data.Box.a22 / cutoff));

    if (cl.Mx < 1)
        cl.Mx = 1;
    if (cl.My < 1)
        cl.My = 1;
    if (cl.Mz < 1)
        cl.Mz = 1;

    cl.Wx = (data.Box.a00 / cl.Mx);
    cl.Wy = (data.Box.a11 / cl.My);
    cl.Wz = (data.Box.a22 / cl.Mz);

    cl.cell_num = cl.Mx * cl.My * cl.Mz;

    cl.Cell.resize(n);
    cl.cell_offset.resize(cl.cell_num + 1);
    cl.atom_indices.resize(n);
}

void lcl1(Data &data, Cell_List &cl, const parameter1 &pr1, const parameter2 &pr2, vector<double> &U_atom)
{
    int n = data.n;
    if (cl.cell_num <= 0 || cl.Mx <= 0 || cl.My <= 0 || cl.Mz <= 0)
    {
        std::cerr << "[lcl1] Invalid cell list.\n";
        return;
    }

    double Wx = cl.Wx;
    double Wy = cl.Wy;
    double Wz = cl.Wz;

    int Mx = cl.Mx;
    int My = cl.My;
    int Mz = cl.Mz;

    auto wrap_index = [](int idx, int m)
    {
        if (m <= 0)
            return 0;
        idx %= m;
        if (idx < 0)
            idx += m;
        return idx;
    };

    std::vector<int> num_in_cell(cl.cell_num, 0);

    for (int i = 0; i < n; i++)
    {
        int mx = static_cast<int>(std::floor(data.atoms[i].r.a00 / Wx));
        int my = static_cast<int>(std::floor(data.atoms[i].r.a10 / Wy));
        int mz = static_cast<int>(std::floor(data.atoms[i].r.a20 / Wz));
        mx = wrap_index(mx, Mx);
        my = wrap_index(my, My);
        mz = wrap_index(mz, Mz);
        cl.Cell[i] = mx + my * Mx + mz * Mx * My; // m is particle index,M is the number of cells
        num_in_cell[cl.Cell[i]]++;
    }

    cl.cell_offset[0] = 0;
    for (int i = 0; i < cl.cell_num; i++)
    {
        cl.cell_offset[i + 1] = cl.cell_offset[i] + num_in_cell[i];
    }

    std::vector<int> offset1(cl.cell_num);
    for (int i = 0; i < cl.cell_num; i++)
    {
        offset1[i] = cl.cell_offset[i];
    }

    for (int i = 0; i < n; i++)
    {
        int c = cl.Cell[i];       // 读取第i个粒子所在的格子编号
        int pos = offset1[c]++;   // 该格子中粒子的个数加1，并返回该粒子在该格子中的位置
        cl.atom_indices[pos] = i; // 将粒子i的编号存入该格子中
    }
    /*
        // 1) 最后一个 offset 必须等于 n
        if (cl.cell_offset[cl.cell_num] != n)
        {
            std::cout << "BAD offset: cell_offset[last]="
                      << cl.cell_offset[cl.cell_num] << " n=" << n << "\n";
            std::abort();
        }

        // 2) offset 必须单调不减，且范围在 [0, n]
        for (int c = 0; c < cl.cell_num; c++)
        {
            if (cl.cell_offset[c] < 0 || cl.cell_offset[c] > n ||
                cl.cell_offset[c + 1] < 0 || cl.cell_offset[c + 1] > n ||
                cl.cell_offset[c] > cl.cell_offset[c + 1])
            {
                std::cout << "BAD offset at c=" << c
                          << " off=" << cl.cell_offset[c]
                          << " off2=" << cl.cell_offset[c + 1] << "\n";
                std::abort();
            }
        }

        // 3) atom_indices 里的内容必须是合法粒子 id
        for (int t = 0; t < n; t++)
        {
            int id = cl.atom_indices[t];
            if (id < 0 || id >= n)
            {
                std::cout << "BAD atom_indices[" << t << "]=" << id << "\n";
                std::abort();
            }
        }
    */
}

void lcl2(Data &data, Cell_List &cl, const parameter1 &pr1,
          const parameter2 &pr2_WW,
          const parameter2 &pr2_BB,
          const parameter2 &pr2_WB,
          vector<double> &U_atom)
{
    int n = data.n;

    int Mx = cl.Mx;
    int My = cl.My;
    int Mz = cl.Mz;
    if (cl.cell_num <= 0 || Mx <= 0 || My <= 0 || Mz <= 0)
    {
        std::cerr << "[lcl2] Invalid cell list.\n";
        return;
    }

    data.F_all = Matrix31(0.0, 0.0, 0.0);
    for (int i = 0; i < n; i++)
    {
        data.atoms[i].f = Matrix31(0.0, 0.0, 0.0);
    }
    U_atom.assign(n, 0.0);

    double Lx = data.Box.a00;
    double Ly = data.Box.a11;
    double Lz = data.Box.a22;
    double Lxh = Lx / 2.0;
    double Lyh = Ly / 2.0;
    double Lzh = Lz / 2.0;

    // Map species strings to small integer tags once per step.
    // 0 = W, 1 = Be, 2 = unknown (fallback to mixed pair parameters).
    std::vector<unsigned char> atom_type(n, 2);
    for (int i = 0; i < n; ++i)
    {
        const std::string &name = data.atoms[i].name;
        atom_type[i] = (name == "W") ? 0 : ((name == "Be") ? 1 : 2);
    }

    // Position cache in SoA layout for hot loops in lcl2.
    std::vector<double> rx_vec(n), ry_vec(n), rz_vec(n);
    for (int i = 0; i < n; ++i)
    {
        rx_vec[i] = data.atoms[i].r.a00;
        ry_vec[i] = data.atoms[i].r.a10;
        rz_vec[i] = data.atoms[i].r.a20;
    }
    const double *rx = rx_vec.data();
    const double *ry = ry_vec.data();
    const double *rz = rz_vec.data();

    struct PairCache
    {
        double D0;
        double R;
        double D;
        double mu;
        double r0;
        double gamma;
        double h_ang;

        double aR;
        double aA;
        double bsR;
        double bsA;

        double R_minus_D2;
        double R_plus_D2;
        double half_pi_over_D;

        double c_ang2;
        double d_ang2;
        double inv_d_ang2;
    };

    auto make_cache = [](const parameter2 &p) -> PairCache
    {
        PairCache c{};
        c.D0 = p.D0;
        c.R = p.R;
        c.D = p.D;
        c.mu = p.mu;
        c.r0 = p.r0;
        c.gamma = p.gamma;
        c.h_ang = p.h;

        const double inv_S_minus_1 = 1.0 / (p.S - 1.0);
        c.aR = p.D0 * inv_S_minus_1;
        c.aA = p.S * p.D0 * inv_S_minus_1;
        c.bsR = p.beta * std::sqrt(2.0 * p.S);
        c.bsA = p.beta * std::sqrt(2.0 / p.S);

        const double R_minus_D = p.R - p.D;
        const double R_plus_D = p.R + p.D;
        c.R_minus_D2 = R_minus_D * R_minus_D;
        c.R_plus_D2 = R_plus_D * R_plus_D;
        c.half_pi_over_D = (p.D != 0.0) ? (0.5 * M_PI / p.D) : 0.0;

        c.c_ang2 = p.c * p.c;
        c.d_ang2 = p.d * p.d;
        c.inv_d_ang2 = (c.d_ang2 != 0.0) ? (1.0 / c.d_ang2) : 0.0;
        return c;
    };

    const PairCache cache_WW = make_cache(pr2_WW);
    const PairCache cache_BB = make_cache(pr2_BB);
    const PairCache cache_WB = make_cache(pr2_WB);

    // ---------------------------------------------
    // tools: cell wrap + iterate 27-cell neighbors
    // ---------------------------------------------
    auto wrap_cell = [](int x, int M)
    {
        if (x < 0)
            return x + M;
        if (x >= M)
            return x - M;
        return x;
    };

    auto for_each_in_27_cells = [&](int center, auto &&func)
    {
        int c0 = cl.Cell[center];

        int cz = c0 / (Mx * My);
        int cc = c0 - cz * (Mx * My);
        int cy = cc / Mx;
        int cx = cc - cy * Mx;

        for (int dcz = -1; dcz <= 1; dcz++)
        {
            for (int dcy = -1; dcy <= 1; dcy++)
            {
                for (int dcx = -1; dcx <= 1; dcx++)
                {
                    int cz2 = wrap_cell(cz + dcz, Mz);
                    int cy2 = wrap_cell(cy + dcy, My);
                    int cx2 = wrap_cell(cx + dcx, Mx);

                    int nc = cx2 + cy2 * Mx + cz2 * Mx * My;

                    int begin = cl.cell_offset[nc];
                    int end = cl.cell_offset[nc + 1];

                    for (int iii = begin; iii < end; iii++)
                    {
                        int pid = cl.atom_indices[iii];
                        func(pid);
                    }
                }
            }
        }
    };

// ---------------------------------------------
// main loops
// ---------------------------------------------
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int max_atoms_per_cell = 0;
    for (int c = 0; c < cl.cell_num; ++c)
    {
        int count = cl.cell_offset[c + 1] - cl.cell_offset[c];
        if (count > max_atoms_per_cell)
            max_atoms_per_cell = count;
    }

    int scratch_reserve = 27 * max_atoms_per_cell;
    if (scratch_reserve < 64)
        scratch_reserve = 64;
    if (scratch_reserve > n)
        scratch_reserve = n;

    struct ThreadScratch
    {
        std::vector<int> k_list_ij;
        std::vector<int> k_list_ji;

        std::vector<double> dXij_dri_x_list, dXij_dri_y_list, dXij_dri_z_list;
        std::vector<double> dXij_drj_x_list, dXij_drj_y_list, dXij_drj_z_list;
        std::vector<double> dXij_drk_x_list, dXij_drk_y_list, dXij_drk_z_list;

        std::vector<double> dXji_dri_x_list, dXji_dri_y_list, dXji_dri_z_list;
        std::vector<double> dXji_drj_x_list, dXji_drj_y_list, dXji_drj_z_list;
        std::vector<double> dXji_drk_x_list, dXji_drk_y_list, dXji_drk_z_list;

        void EnsureCapacity(size_t cap)
        {
            auto ensure = [cap](auto &v)
            {
                if (v.capacity() < cap)
                    v.reserve(cap);
            };
            ensure(k_list_ij);
            ensure(k_list_ji);
            ensure(dXij_dri_x_list);
            ensure(dXij_dri_y_list);
            ensure(dXij_dri_z_list);
            ensure(dXij_drj_x_list);
            ensure(dXij_drj_y_list);
            ensure(dXij_drj_z_list);
            ensure(dXij_drk_x_list);
            ensure(dXij_drk_y_list);
            ensure(dXij_drk_z_list);
            ensure(dXji_dri_x_list);
            ensure(dXji_dri_y_list);
            ensure(dXji_dri_z_list);
            ensure(dXji_drj_x_list);
            ensure(dXji_drj_y_list);
            ensure(dXji_drj_z_list);
            ensure(dXji_drk_x_list);
            ensure(dXji_drk_y_list);
            ensure(dXji_drk_z_list);
        }
    };

    std::vector<double> fx_vec(n, 0.0);
    std::vector<double> fy_vec(n, 0.0);
    std::vector<double> fz_vec(n, 0.0);
    std::vector<double> u_vec(n, 0.0);
    double *fx = fx_vec.data();
    double *fy = fy_vec.data();
    double *fz = fz_vec.data();
    double *u = u_vec.data();

#pragma omp parallel
    {
        thread_local ThreadScratch tls;
        tls.EnsureCapacity(static_cast<size_t>(scratch_reserve));

        auto &k_list_ij = tls.k_list_ij;
        auto &k_list_ji = tls.k_list_ji;

        auto &dXij_dri_x_list = tls.dXij_dri_x_list;
        auto &dXij_dri_y_list = tls.dXij_dri_y_list;
        auto &dXij_dri_z_list = tls.dXij_dri_z_list;
        auto &dXij_drj_x_list = tls.dXij_drj_x_list;
        auto &dXij_drj_y_list = tls.dXij_drj_y_list;
        auto &dXij_drj_z_list = tls.dXij_drj_z_list;
        auto &dXij_drk_x_list = tls.dXij_drk_x_list;
        auto &dXij_drk_y_list = tls.dXij_drk_y_list;
        auto &dXij_drk_z_list = tls.dXij_drk_z_list;

        auto &dXji_dri_x_list = tls.dXji_dri_x_list;
        auto &dXji_dri_y_list = tls.dXji_dri_y_list;
        auto &dXji_dri_z_list = tls.dXji_dri_z_list;
        auto &dXji_drj_x_list = tls.dXji_drj_x_list;
        auto &dXji_drj_y_list = tls.dXji_drj_y_list;
        auto &dXji_drj_z_list = tls.dXji_drj_z_list;
        auto &dXji_drk_x_list = tls.dXji_drk_x_list;
        auto &dXji_drk_y_list = tls.dXji_drk_y_list;
        auto &dXji_drk_z_list = tls.dXji_drk_z_list;

#pragma omp for schedule(guided, 32) reduction(+:fx[:n], fy[:n], fz[:n], u[:n])
        for (int i = 0; i < n; i++)
        {
            const unsigned char type_i = atom_type[i];

            int c_i = cl.Cell[i];
            int cz_i = c_i / (Mx * My);
            int cc_i = c_i - cz_i * (Mx * My); // just a number
            int cy_i = cc_i / Mx;
            int cx_i = cc_i - cy_i * Mx;

            auto process_pair = [&](int j)
            {
                const unsigned char type_j = atom_type[j];

                const PairCache *cache = &cache_WB;
                if (type_i == 0 && type_j == 0)
                {
                    cache = &cache_WW;
                }
                else if (type_i == 1 && type_j == 1)
                {
                    cache = &cache_BB;
                }

                const double mu = cache->mu;
                const double r0 = cache->r0;
                const double gamma = cache->gamma;
                const double h_ang = cache->h_ang;

                // displacement r_ij = r_j - r_i (scalar SoA-friendly path)
                double drij_x = rx[j] - rx[i];
                double drij_y = ry[j] - ry[i];
                double drij_z = rz[j] - rz[i];

                // PBC (real space)
                if (drij_x < -Lxh)
                    drij_x += Lx;
                else if (drij_x >= Lxh)
                    drij_x -= Lx;

                if (drij_y < -Lyh)
                    drij_y += Ly;
                else if (drij_y >= Lyh)
                    drij_y -= Ly;

                if (drij_z < -Lzh)
                    drij_z += Lz;
                else if (drij_z >= Lzh)
                    drij_z -= Lz;

                double rij2 = drij_x * drij_x + drij_y * drij_y + drij_z * drij_z;
                if (rij2 == 0.0)
                {
                    return;
                }
                if (rij2 >= cache->R_plus_D2)
                {
                    return;
                }

                double rij = std::sqrt(rij2);
                double inv_rij = 1.0 / rij;
                double seij_x = drij_x * inv_rij;
                double seij_y = drij_y * inv_rij;
                double seij_z = drij_z * inv_rij;
                double drji_x = -drij_x;
                double drji_y = -drij_y;
                double drji_z = -drij_z;
                double rji = rij;
                double seji_x = -seij_x;
                double seji_y = -seij_y;
                double seji_z = -seij_z;

                // Pair VR / VA
                const double expR = std::exp(-cache->bsR * (rij - r0));
                const double expA = std::exp(-cache->bsA * (rij - r0));
                const double VR = cache->aR * expR;
                const double VA = cache->aA * expA;

                // Cutoff fc_ij and derivative
                double fcij = 0.0;
                double nfcij = 0.0;
                if (rij2 <= cache->R_minus_D2)
                {
                    fcij = 1.0;
                    nfcij = 0.0;
                }
                else
                {
                    const double x = (rij - cache->R) * cache->half_pi_over_D;
                    fcij = 0.5 - 0.5 * std::sin(x);
                    nfcij = -0.5 * std::cos(x) * cache->half_pi_over_D;
                }

                if (fcij == 0.0)
                {
                    return;
                }

                const double nVR = -cache->bsR * VR;
                const double nVA = -cache->bsA * VA;

                // -------------------------------------------------
                // Three-body: build Xij using k in i-neighborhood,
                //            build Xji using k in j-neighborhood
                // -------------------------------------------------
                double Xij = 0.0;
                double Xji = 0.0;

                k_list_ij.clear();
                k_list_ji.clear();

                dXij_dri_x_list.clear();
                dXij_dri_y_list.clear();
                dXij_dri_z_list.clear();
                dXij_drj_x_list.clear();
                dXij_drj_y_list.clear();
                dXij_drj_z_list.clear();
                dXij_drk_x_list.clear();
                dXij_drk_y_list.clear();
                dXij_drk_z_list.clear();

                dXji_dri_x_list.clear();
                dXji_dri_y_list.clear();
                dXji_dri_z_list.clear();
                dXji_drj_x_list.clear();
                dXji_drj_y_list.clear();
                dXji_drj_z_list.clear();
                dXji_drk_x_list.clear();
                dXji_drk_y_list.clear();
                dXji_drk_z_list.clear();

                auto accumulate_X_center = [&](int center,
                                               double ac_x, double ac_y, double ac_z,
                                               double dr_co_x, double dr_co_y, double dr_co_z,
                                               double se_co_x, double se_co_y, double se_co_z,
                                               double r_co,
                                               double &X,
                                               std::vector<double> &dX_dr_center_x_list,
                                               std::vector<double> &dX_dr_center_y_list,
                                               std::vector<double> &dX_dr_center_z_list,
                                               std::vector<double> &dX_dr_other_x_list,
                                               std::vector<double> &dX_dr_other_y_list,
                                               std::vector<double> &dX_dr_other_z_list,
                                               std::vector<double> &dX_drk_x_list,
                                               std::vector<double> &dX_drk_y_list,
                                               std::vector<double> &dX_drk_z_list,
                                               std::vector<int> &k_list)
                {
                    const double inv_r_co = 1.0 / r_co;
                    const double inv_r_co2 = inv_r_co * inv_r_co;

                    for_each_in_27_cells(center, [&](int k)
                    {
                        if (k == i || k == j)
                        {
                            return;
                        }

                        double dr_ck_x = rx[k] - ac_x;
                        double dr_ck_y = ry[k] - ac_y;
                        double dr_ck_z = rz[k] - ac_z;

                        if (dr_ck_x < -Lxh)
                            dr_ck_x += Lx;
                        else if (dr_ck_x >= Lxh)
                            dr_ck_x -= Lx;

                        if (dr_ck_y < -Lyh)
                            dr_ck_y += Ly;
                        else if (dr_ck_y >= Lyh)
                            dr_ck_y -= Ly;

                        if (dr_ck_z < -Lzh)
                            dr_ck_z += Lz;
                        else if (dr_ck_z >= Lzh)
                            dr_ck_z -= Lz;

                        double r_ck2 = dr_ck_x * dr_ck_x + dr_ck_y * dr_ck_y + dr_ck_z * dr_ck_z;
                        if (r_ck2 == 0.0)
                        {
                            return;
                        }
                        if (r_ck2 >= cache->R_plus_D2)
                        {
                            return;
                        }
                        double r_ck = std::sqrt(r_ck2);
                        double inv_r_ck = 1.0 / r_ck;
                        double inv_r_ck2 = inv_r_ck * inv_r_ck;
                        double inv_r_co_ck = inv_r_co * inv_r_ck;

                        double fc_ck, dfc_dr_ck;
                        if (r_ck2 <= cache->R_minus_D2)
                        {
                            fc_ck = 1.0;
                            dfc_dr_ck = 0.0;
                        }
                        else
                        {
                            const double xx = (r_ck - cache->R) * cache->half_pi_over_D;
                            fc_ck = 0.5 - 0.5 * std::sin(xx);
                            dfc_dr_ck = -0.5 * std::cos(xx) * cache->half_pi_over_D;
                        }

                        if (fc_ck <= 0.0)
                        {
                            return;
                        }

                        double se_ck_x = dr_ck_x * inv_r_ck;
                        double se_ck_y = dr_ck_y * inv_r_ck;
                        double se_ck_z = dr_ck_z * inv_r_ck;

                        double dot_co_ck = dr_co_x * dr_ck_x + dr_co_y * dr_ck_y + dr_co_z * dr_ck_z;
                        double cost_cok = dot_co_ck * inv_r_co_ck;
                        if (cost_cok > 1.0)
                            cost_cok = 1.0;
                        if (cost_cok < -1.0)
                            cost_cok = -1.0;

                        const double hc = h_ang + cost_cok;
                        const double u_ang = cache->d_ang2 + hc * hc;
                        const double inv_u_ang = 1.0 / u_ang;
                        const double g_ck = gamma * (1.0 + cache->c_ang2 * (cache->inv_d_ang2 - inv_u_ang));
                        const double dg_dcos = gamma * cache->c_ang2 * 2.0 * hc * inv_u_ang * inv_u_ang;

                        double term_dr_co = dot_co_ck * inv_r_co2;
                        double dcos_ddr_co_x = (dr_ck_x - dr_co_x * term_dr_co) * inv_r_co_ck;
                        double dcos_ddr_co_y = (dr_ck_y - dr_co_y * term_dr_co) * inv_r_co_ck;
                        double dcos_ddr_co_z = (dr_ck_z - dr_co_z * term_dr_co) * inv_r_co_ck;

                        double term_dr_ck = dot_co_ck * inv_r_ck2;
                        double dcos_ddr_ck_x = (dr_co_x - dr_ck_x * term_dr_ck) * inv_r_co_ck;
                        double dcos_ddr_ck_y = (dr_co_y - dr_ck_y * term_dr_ck) * inv_r_co_ck;
                        double dcos_ddr_ck_z = (dr_co_z - dr_ck_z * term_dr_ck) * inv_r_co_ck;

                        double dcos_drc_x = -(dcos_ddr_co_x + dcos_ddr_ck_x);
                        double dcos_drc_y = -(dcos_ddr_co_y + dcos_ddr_ck_y);
                        double dcos_drc_z = -(dcos_ddr_co_z + dcos_ddr_ck_z);
                        double dcos_dro_x = dcos_ddr_co_x;
                        double dcos_dro_y = dcos_ddr_co_y;
                        double dcos_dro_z = dcos_ddr_co_z;
                        double dcos_drk_x = dcos_ddr_ck_x;
                        double dcos_drk_y = dcos_ddr_ck_y;
                        double dcos_drk_z = dcos_ddr_ck_z;

                        double dg_drc_x = dcos_drc_x * dg_dcos;
                        double dg_drc_y = dcos_drc_y * dg_dcos;
                        double dg_drc_z = dcos_drc_z * dg_dcos;
                        double dg_dro_x = dcos_dro_x * dg_dcos;
                        double dg_dro_y = dcos_dro_y * dg_dcos;
                        double dg_dro_z = dcos_dro_z * dg_dcos;
                        double dg_drk_x = dcos_drk_x * dg_dcos;
                        double dg_drk_y = dcos_drk_y * dg_dcos;
                        double dg_drk_z = dcos_drk_z * dg_dcos;

                        double e_cok = std::exp(mu * (r_co - r_ck));
                        double mu_e = mu * e_cok;

                        double dfc_drc_x = -se_ck_x * dfc_dr_ck;
                        double dfc_drc_y = -se_ck_y * dfc_dr_ck;
                        double dfc_drc_z = -se_ck_z * dfc_dr_ck;
                        double dfc_drk_x = se_ck_x * dfc_dr_ck;
                        double dfc_drk_y = se_ck_y * dfc_dr_ck;
                        double dfc_drk_z = se_ck_z * dfc_dr_ck;

                        double de_drc_x = (se_ck_x - se_co_x) * mu_e;
                        double de_drc_y = (se_ck_y - se_co_y) * mu_e;
                        double de_drc_z = (se_ck_z - se_co_z) * mu_e;
                        double de_dro_x = se_co_x * mu_e;
                        double de_dro_y = se_co_y * mu_e;
                        double de_dro_z = se_co_z * mu_e;
                        double de_drk_x = -se_ck_x * mu_e;
                        double de_drk_y = -se_ck_y * mu_e;
                        double de_drk_z = -se_ck_z * mu_e;

                        double g_e = g_ck * e_cok;
                        double fc_e = fc_ck * e_cok;
                        double fc_g = fc_ck * g_ck;

                        double dX_drc_x = dfc_drc_x * g_e + dg_drc_x * fc_e + de_drc_x * fc_g;
                        double dX_drc_y = dfc_drc_y * g_e + dg_drc_y * fc_e + de_drc_y * fc_g;
                        double dX_drc_z = dfc_drc_z * g_e + dg_drc_z * fc_e + de_drc_z * fc_g;
                        double dX_dro_x = dg_dro_x * fc_e + de_dro_x * fc_g;
                        double dX_dro_y = dg_dro_y * fc_e + de_dro_y * fc_g;
                        double dX_dro_z = dg_dro_z * fc_e + de_dro_z * fc_g;
                        double dX_drk_x = dfc_drk_x * g_e + dg_drk_x * fc_e + de_drk_x * fc_g;
                        double dX_drk_y = dfc_drk_y * g_e + dg_drk_y * fc_e + de_drk_y * fc_g;
                        double dX_drk_z = dfc_drk_z * g_e + dg_drk_z * fc_e + de_drk_z * fc_g;

                        X += fc_ck * g_ck * e_cok;

                        dX_dr_center_x_list.push_back(dX_drc_x);
                        dX_dr_center_y_list.push_back(dX_drc_y);
                        dX_dr_center_z_list.push_back(dX_drc_z);
                        dX_dr_other_x_list.push_back(dX_dro_x);
                        dX_dr_other_y_list.push_back(dX_dro_y);
                        dX_dr_other_z_list.push_back(dX_dro_z);
                        dX_drk_x_list.push_back(dX_drk_x);
                        dX_drk_y_list.push_back(dX_drk_y);
                        dX_drk_z_list.push_back(dX_drk_z);
                        k_list.push_back(k);
                    });
                };

                // Xij : i-centered
                accumulate_X_center(i, rx[i], ry[i], rz[i],
                                    drij_x, drij_y, drij_z,
                                    seij_x, seij_y, seij_z,
                                    rij, Xij,
                                    dXij_dri_x_list, dXij_dri_y_list, dXij_dri_z_list,
                                    dXij_drj_x_list, dXij_drj_y_list, dXij_drj_z_list,
                                    dXij_drk_x_list, dXij_drk_y_list, dXij_drk_z_list,
                                    k_list_ij);

                // Xji : j-centered
                accumulate_X_center(j, rx[j], ry[j], rz[j],
                                    drji_x, drji_y, drji_z,
                                    seji_x, seji_y, seji_z,
                                    rji, Xji,
                                    dXji_drj_x_list, dXji_drj_y_list, dXji_drj_z_list,
                                    dXji_dri_x_list, dXji_dri_y_list, dXji_dri_z_list,
                                    dXji_drk_x_list, dXji_drk_y_list, dXji_drk_z_list,
                                    k_list_ji);

                // bond order terms
                double bij = 1.0 / std::sqrt(1.0 + Xij);
                double bji = 1.0 / std::sqrt(1.0 + Xji);
                double b_sym = 0.5 * (bij + bji);

                double u_pair = fcij * (VR - b_sym * VA);
                u[i] += 0.5 * u_pair;
                u[j] += 0.5 * u_pair;

                double denom_ij = std::sqrt(1.0 + Xij) * (1.0 + Xij);
                double denom_ji = std::sqrt(1.0 + Xji) * (1.0 + Xji);

                double fk_ij = 0.0;
                double fk_ji = 0.0;

                if (denom_ij > 0.0)
                {
                    fk_ij = fcij * VA * (0.25 / denom_ij);
                }
                if (denom_ji > 0.0)
                {
                    fk_ji = fcij * VA * (0.25 / denom_ji);
                }

                // forces from Xij (i-centered)
                for (size_t idx = 0; idx < k_list_ij.size(); ++idx)
                {
                    int kk = k_list_ij[idx];
                    fx[i] += -fk_ij * dXij_dri_x_list[idx];
                    fy[i] += -fk_ij * dXij_dri_y_list[idx];
                    fz[i] += -fk_ij * dXij_dri_z_list[idx];

                    fx[j] += -fk_ij * dXij_drj_x_list[idx];
                    fy[j] += -fk_ij * dXij_drj_y_list[idx];
                    fz[j] += -fk_ij * dXij_drj_z_list[idx];

                    fx[kk] += -fk_ij * dXij_drk_x_list[idx];
                    fy[kk] += -fk_ij * dXij_drk_y_list[idx];
                    fz[kk] += -fk_ij * dXij_drk_z_list[idx];
                }

                // forces from Xji (j-centered)
                for (size_t idx = 0; idx < k_list_ji.size(); ++idx)
                {
                    int kk = k_list_ji[idx];
                    fx[i] += -fk_ji * dXji_dri_x_list[idx];
                    fy[i] += -fk_ji * dXji_dri_y_list[idx];
                    fz[i] += -fk_ji * dXji_dri_z_list[idx];

                    fx[j] += -fk_ji * dXji_drj_x_list[idx];
                    fy[j] += -fk_ji * dXji_drj_y_list[idx];
                    fz[j] += -fk_ji * dXji_drj_z_list[idx];

                    fx[kk] += -fk_ji * dXji_drk_x_list[idx];
                    fy[kk] += -fk_ji * dXji_drk_y_list[idx];
                    fz[kk] += -fk_ji * dXji_drk_z_list[idx];
                }

                // pair (radial) forces from fc, VR, VA
                double Fij_scalar = nfcij * (VR - b_sym * VA) + fcij * (nVR - b_sym * nVA);

                fx[i] += seij_x * Fij_scalar;
                fy[i] += seij_y * Fij_scalar;
                fz[i] += seij_z * Fij_scalar;
                fx[j] -= seij_x * Fij_scalar;
                fy[j] -= seij_y * Fij_scalar;
                fz[j] -= seij_z * Fij_scalar;
            };

            // Toggle traversal mode by commenting one of the two lines below.
            // const bool use_half_stencil = true;
            const bool use_half_stencil = false;//!(14 ture <-> false 27)

            if (use_half_stencil)
            {
                static const int half_offsets[14][3] = {
                    {0, 0, 0},
                    {1, 0, 0},
                    {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
                    {-1, -1, 1}, {0, -1, 1}, {1, -1, 1},
                    {-1, 0, 1}, {0, 0, 1}, {1, 0, 1},
                    {-1, 1, 1}, {0, 1, 1}, {1, 1, 1}};

                int neighbor_cells[14];
                int neighbor_count = 0;
                for (int h = 0; h < 14; ++h)
                {
                    int cz2 = wrap_cell(cz_i + half_offsets[h][2], Mz);
                    int cy2 = wrap_cell(cy_i + half_offsets[h][1], My);
                    int cx2 = wrap_cell(cx_i + half_offsets[h][0], Mx);
                    int nc = cx2 + cy2 * Mx + cz2 * Mx * My;

                    bool duplicate = false;
                    for (int t = 0; t < neighbor_count; ++t)
                    {
                        if (neighbor_cells[t] == nc)
                        {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate)
                    {
                        neighbor_cells[neighbor_count++] = nc;
                    }
                }

                for (int h = 0; h < neighbor_count; ++h)
                {
                    int nc = neighbor_cells[h];
                    const bool same_cell = (nc == c_i);

                    int begin = cl.cell_offset[nc];
                    int end = cl.cell_offset[nc + 1];
                    for (int jj = begin; jj < end; ++jj)
                    {
                        int j = cl.atom_indices[jj];
                        if (same_cell && j <= i)
                        {
                            continue;
                        }
                        process_pair(j);
                    }
                }
            }
            else
            {
                for (int dcz = -1; dcz <= 1; ++dcz)
                {
                    for (int dcy = -1; dcy <= 1; ++dcy)
                    {
                        for (int dcx = -1; dcx <= 1; ++dcx)
                        {
                            int cz2 = wrap_cell(cz_i + dcz, Mz);
                            int cy2 = wrap_cell(cy_i + dcy, My);
                            int cx2 = wrap_cell(cx_i + dcx, Mx);
                            int nc = cx2 + cy2 * Mx + cz2 * Mx * My;

                            int begin = cl.cell_offset[nc];
                            int end = cl.cell_offset[nc + 1];
                            for (int jj = begin; jj < end; ++jj)
                            {
                                int j = cl.atom_indices[jj];
                                if (j <= i)
                                {
                                    continue;
                                }
                                process_pair(j);
                            }
                        }
                    }
                }
            }
        }

    }

    for (int i = 0; i < n; i++)
    {
        data.atoms[i].f.a00 = fx[i];
        data.atoms[i].f.a10 = fy[i];
        data.atoms[i].f.a20 = fz[i];
        U_atom[i] = u[i];

        data.F_all.a00 = data.F_all.a00 + data.atoms[i].f.a00;
        data.F_all.a10 = data.F_all.a10 + data.atoms[i].f.a10;
        data.F_all.a20 = data.F_all.a20 + data.atoms[i].f.a20;
    }
}
