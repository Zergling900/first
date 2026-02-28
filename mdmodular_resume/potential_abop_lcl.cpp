#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

#include "3.h"
#include "void.h"

void lcl2(Data &data, Cell_List &cl, const parameter1 &pr1,
          const parameter2 &pr2_WW,
          const parameter2 &pr2_BB,
          const parameter2 &pr2_WB,
          vector<double> &U_atom)
{
    (void)pr1;
    const int n = data.n;

    if (!cl.verlet_valid || static_cast<int>(cl.verlet_offset.size()) != n + 1)
    {
        std::cerr << "[lcl2] Invalid Verlet list.\n";
        return;
    }

    data.F_all = Matrix31(0.0, 0.0, 0.0);
    U_atom.assign(n, 0.0);

    const double Lx = data.Box.a00;
    const double Ly = data.Box.a11;
    const double Lz = data.Box.a22;
    const double Lxh = 0.5 * Lx;
    const double Lyh = 0.5 * Ly;
    const double Lzh = 0.5 * Lz;

    // SoA cache for hot loop access.
    std::vector<double> rx_vec(n), ry_vec(n), rz_vec(n);
    std::vector<unsigned char> type_vec(n, ATOM_TYPE_OTHER);
    for (int i = 0; i < n; ++i)
    {
        rx_vec[i] = data.atoms[i].r.a00;
        ry_vec[i] = data.atoms[i].r.a10;
        rz_vec[i] = data.atoms[i].r.a20;
        const unsigned char t = data.atoms[i].atom_type;
        type_vec[i] = (t <= ATOM_TYPE_OTHER) ? t : ATOM_TYPE_OTHER;
    }
    const double *rx = rx_vec.data();
    const double *ry = ry_vec.data();
    const double *rz = rz_vec.data();
    const unsigned char *types = type_vec.data();

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

    // Branchless pair cache dispatch table.
    const PairCache *cache_lut[3][3] = {
        {&cache_WW, &cache_WB, &cache_WB},
        {&cache_WB, &cache_BB, &cache_WB},
        {&cache_WB, &cache_WB, &cache_WB}};

    int max_neighbors = 0;
    for (int i = 0; i < n; ++i)
    {
        const int count = cl.verlet_offset[i + 1] - cl.verlet_offset[i];
        if (count > max_neighbors)
            max_neighbors = count;
    }
    int scratch_size = max_neighbors;
    if (scratch_size < 64)
        scratch_size = 64;
    if (scratch_size > n)
        scratch_size = n;

    struct DerivBuf
    {
        std::vector<int> k_list;
        std::vector<double> d_center_x, d_center_y, d_center_z;
        std::vector<double> d_other_x, d_other_y, d_other_z;
        std::vector<double> d_k_x, d_k_y, d_k_z;

        void EnsureSize(size_t cap)
        {
            auto ensure = [cap](auto &v)
            {
                if (v.size() < cap)
                    v.resize(cap);
            };
            ensure(k_list);
            ensure(d_center_x);
            ensure(d_center_y);
            ensure(d_center_z);
            ensure(d_other_x);
            ensure(d_other_y);
            ensure(d_other_z);
            ensure(d_k_x);
            ensure(d_k_y);
            ensure(d_k_z);
        }
    };

    struct ThreadScratch
    {
        DerivBuf ij;
        DerivBuf ji;

        void EnsureSize(size_t cap)
        {
            ij.EnsureSize(cap);
            ji.EnsureSize(cap);
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
        tls.EnsureSize(static_cast<size_t>(scratch_size));

        auto ensure_slot = [&](DerivBuf &buf, int idx)
        {
            if (idx < static_cast<int>(buf.k_list.size()))
                return;
            size_t new_cap = buf.k_list.empty() ? static_cast<size_t>(64) : buf.k_list.size();
            while (idx >= static_cast<int>(new_cap))
                new_cap = new_cap * 2;
            buf.EnsureSize(new_cap);
        };

#pragma omp for schedule(guided, 32) reduction(+ : fx[:n], fy[:n], fz[:n], u[:n])
        for (int i = 0; i < n; ++i)
        {
            const unsigned char type_i = types[i];

            const int begin = cl.verlet_offset[i];
            const int end = cl.verlet_offset[i + 1];
            for (int jj = begin; jj < end; ++jj)
            {
                const int j = cl.verlet_neighbors[jj];
                if (j <= i)
                    continue;

                const unsigned char type_j = types[j];
                const PairCache &cache = *cache_lut[type_i][type_j];

                const double mu = cache.mu;
                const double r0 = cache.r0;
                const double gamma = cache.gamma;
                const double h_ang = cache.h_ang;

                double drij_x = rx[j] - rx[i];
                double drij_y = ry[j] - ry[i];
                double drij_z = rz[j] - rz[i];

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

                const double rij2 = drij_x * drij_x + drij_y * drij_y + drij_z * drij_z;
                if (rij2 == 0.0 || rij2 >= cache.R_plus_D2)
                    continue;

                const double rij = std::sqrt(rij2);
                const double inv_rij = 1.0 / rij;
                const double seij_x = drij_x * inv_rij;
                const double seij_y = drij_y * inv_rij;
                const double seij_z = drij_z * inv_rij;
                const double drji_x = -drij_x;
                const double drji_y = -drij_y;
                const double drji_z = -drij_z;
                const double seji_x = -seij_x;
                const double seji_y = -seij_y;
                const double seji_z = -seij_z;

                const double expR = std::exp(-cache.bsR * (rij - r0));
                const double expA = std::exp(-cache.bsA * (rij - r0));
                const double VR = cache.aR * expR;
                const double VA = cache.aA * expA;

                double fcij = 0.0;
                double nfcij = 0.0;
                if (rij2 <= cache.R_minus_D2)
                {
                    fcij = 1.0;
                }
                else
                {
                    const double x = (rij - cache.R) * cache.half_pi_over_D;
                    fcij = 0.5 - 0.5 * std::sin(x);
                    nfcij = -0.5 * std::cos(x) * cache.half_pi_over_D;
                }
                if (fcij == 0.0)
                    continue;

                const double nVR = -cache.bsR * VR;
                const double nVA = -cache.bsA * VA;

                auto accumulate_X_center = [&](int center,
                                               double ac_x, double ac_y, double ac_z,
                                               double dr_co_x, double dr_co_y, double dr_co_z,
                                               double se_co_x, double se_co_y, double se_co_z,
                                               double r_co,
                                               DerivBuf &buf,
                                               int &term_count,
                                               double &X)
                {
                    const double inv_r_co = 1.0 / r_co;
                    const double inv_r_co2 = inv_r_co * inv_r_co;

                    const int begin_n = cl.verlet_offset[center];
                    const int end_n = cl.verlet_offset[center + 1];
                    for (int nn = begin_n; nn < end_n; ++nn)
                    {
                        const int k = cl.verlet_neighbors[nn];
                        if (k == i || k == j)
                            continue;

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

                        const double r_ck2 = dr_ck_x * dr_ck_x + dr_ck_y * dr_ck_y + dr_ck_z * dr_ck_z;
                        if (r_ck2 == 0.0 || r_ck2 >= cache.R_plus_D2)
                            continue;

                        const double r_ck = std::sqrt(r_ck2);
                        const double inv_r_ck = 1.0 / r_ck;
                        const double inv_r_ck2 = inv_r_ck * inv_r_ck;
                        const double inv_r_co_ck = inv_r_co * inv_r_ck;

                        double fc_ck = 0.0;
                        double dfc_dr_ck = 0.0;
                        if (r_ck2 <= cache.R_minus_D2)
                        {
                            fc_ck = 1.0;
                        }
                        else
                        {
                            const double xx = (r_ck - cache.R) * cache.half_pi_over_D;
                            fc_ck = 0.5 - 0.5 * std::sin(xx);
                            dfc_dr_ck = -0.5 * std::cos(xx) * cache.half_pi_over_D;
                        }
                        if (fc_ck <= 0.0)
                            continue;

                        const double se_ck_x = dr_ck_x * inv_r_ck;
                        const double se_ck_y = dr_ck_y * inv_r_ck;
                        const double se_ck_z = dr_ck_z * inv_r_ck;

                        const double dot_co_ck = dr_co_x * dr_ck_x + dr_co_y * dr_ck_y + dr_co_z * dr_ck_z;
                        double cost_cok = dot_co_ck * inv_r_co_ck;
                        if (cost_cok > 1.0)
                            cost_cok = 1.0;
                        if (cost_cok < -1.0)
                            cost_cok = -1.0;

                        const double hc = h_ang + cost_cok;
                        const double u_ang = cache.d_ang2 + hc * hc;
                        const double inv_u_ang = 1.0 / u_ang;
                        const double g_ck = gamma * (1.0 + cache.c_ang2 * (cache.inv_d_ang2 - inv_u_ang));
                        const double dg_dcos = gamma * cache.c_ang2 * 2.0 * hc * inv_u_ang * inv_u_ang;

                        const double term_dr_co = dot_co_ck * inv_r_co2;
                        const double dcos_ddr_co_x = (dr_ck_x - dr_co_x * term_dr_co) * inv_r_co_ck;
                        const double dcos_ddr_co_y = (dr_ck_y - dr_co_y * term_dr_co) * inv_r_co_ck;
                        const double dcos_ddr_co_z = (dr_ck_z - dr_co_z * term_dr_co) * inv_r_co_ck;

                        const double term_dr_ck = dot_co_ck * inv_r_ck2;
                        const double dcos_ddr_ck_x = (dr_co_x - dr_ck_x * term_dr_ck) * inv_r_co_ck;
                        const double dcos_ddr_ck_y = (dr_co_y - dr_ck_y * term_dr_ck) * inv_r_co_ck;
                        const double dcos_ddr_ck_z = (dr_co_z - dr_ck_z * term_dr_ck) * inv_r_co_ck;

                        const double dcos_drc_x = -(dcos_ddr_co_x + dcos_ddr_ck_x);
                        const double dcos_drc_y = -(dcos_ddr_co_y + dcos_ddr_ck_y);
                        const double dcos_drc_z = -(dcos_ddr_co_z + dcos_ddr_ck_z);
                        const double dcos_dro_x = dcos_ddr_co_x;
                        const double dcos_dro_y = dcos_ddr_co_y;
                        const double dcos_dro_z = dcos_ddr_co_z;
                        const double dcos_drk_x = dcos_ddr_ck_x;
                        const double dcos_drk_y = dcos_ddr_ck_y;
                        const double dcos_drk_z = dcos_ddr_ck_z;

                        const double dg_drc_x = dcos_drc_x * dg_dcos;
                        const double dg_drc_y = dcos_drc_y * dg_dcos;
                        const double dg_drc_z = dcos_drc_z * dg_dcos;
                        const double dg_dro_x = dcos_dro_x * dg_dcos;
                        const double dg_dro_y = dcos_dro_y * dg_dcos;
                        const double dg_dro_z = dcos_dro_z * dg_dcos;
                        const double dg_drk_x = dcos_drk_x * dg_dcos;
                        const double dg_drk_y = dcos_drk_y * dg_dcos;
                        const double dg_drk_z = dcos_drk_z * dg_dcos;

                        const double e_cok = std::exp(mu * (r_co - r_ck));
                        const double mu_e = mu * e_cok;

                        const double dfc_drc_x = -se_ck_x * dfc_dr_ck;
                        const double dfc_drc_y = -se_ck_y * dfc_dr_ck;
                        const double dfc_drc_z = -se_ck_z * dfc_dr_ck;
                        const double dfc_drk_x = se_ck_x * dfc_dr_ck;
                        const double dfc_drk_y = se_ck_y * dfc_dr_ck;
                        const double dfc_drk_z = se_ck_z * dfc_dr_ck;

                        const double de_drc_x = (se_ck_x - se_co_x) * mu_e;
                        const double de_drc_y = (se_ck_y - se_co_y) * mu_e;
                        const double de_drc_z = (se_ck_z - se_co_z) * mu_e;
                        const double de_dro_x = se_co_x * mu_e;
                        const double de_dro_y = se_co_y * mu_e;
                        const double de_dro_z = se_co_z * mu_e;
                        const double de_drk_x = -se_ck_x * mu_e;
                        const double de_drk_y = -se_ck_y * mu_e;
                        const double de_drk_z = -se_ck_z * mu_e;

                        const double g_e = g_ck * e_cok;
                        const double fc_e = fc_ck * e_cok;
                        const double fc_g = fc_ck * g_ck;

                        const double dX_drc_x = dfc_drc_x * g_e + dg_drc_x * fc_e + de_drc_x * fc_g;
                        const double dX_drc_y = dfc_drc_y * g_e + dg_drc_y * fc_e + de_drc_y * fc_g;
                        const double dX_drc_z = dfc_drc_z * g_e + dg_drc_z * fc_e + de_drc_z * fc_g;
                        const double dX_dro_x = dg_dro_x * fc_e + de_dro_x * fc_g;
                        const double dX_dro_y = dg_dro_y * fc_e + de_dro_y * fc_g;
                        const double dX_dro_z = dg_dro_z * fc_e + de_dro_z * fc_g;
                        const double dX_drk_x = dfc_drk_x * g_e + dg_drk_x * fc_e + de_drk_x * fc_g;
                        const double dX_drk_y = dfc_drk_y * g_e + dg_drk_y * fc_e + de_drk_y * fc_g;
                        const double dX_drk_z = dfc_drk_z * g_e + dg_drk_z * fc_e + de_drk_z * fc_g;

                        X += fc_ck * g_ck * e_cok;

                        ensure_slot(buf, term_count);
                        buf.k_list[term_count] = k;
                        buf.d_center_x[term_count] = dX_drc_x;
                        buf.d_center_y[term_count] = dX_drc_y;
                        buf.d_center_z[term_count] = dX_drc_z;
                        buf.d_other_x[term_count] = dX_dro_x;
                        buf.d_other_y[term_count] = dX_dro_y;
                        buf.d_other_z[term_count] = dX_dro_z;
                        buf.d_k_x[term_count] = dX_drk_x;
                        buf.d_k_y[term_count] = dX_drk_y;
                        buf.d_k_z[term_count] = dX_drk_z;
                        ++term_count;
                    }
                };

                DerivBuf &buf_ij = tls.ij;
                DerivBuf &buf_ji = tls.ji;
                int count_ij = 0;
                int count_ji = 0;
                double Xij = 0.0;
                double Xji = 0.0;

                accumulate_X_center(i,
                                    rx[i], ry[i], rz[i],
                                    drij_x, drij_y, drij_z,
                                    seij_x, seij_y, seij_z,
                                    rij,
                                    buf_ij,
                                    count_ij,
                                    Xij);

                accumulate_X_center(j,
                                    rx[j], ry[j], rz[j],
                                    drji_x, drji_y, drji_z,
                                    seji_x, seji_y, seji_z,
                                    rij,
                                    buf_ji,
                                    count_ji,
                                    Xji);

                const double one_p_Xij = 1.0 + Xij;
                const double one_p_Xji = 1.0 + Xji;
                const double sqrt_one_p_Xij = std::sqrt(one_p_Xij);
                const double sqrt_one_p_Xji = std::sqrt(one_p_Xji);

                const double bij = 1.0 / sqrt_one_p_Xij;
                const double bji = 1.0 / sqrt_one_p_Xji;
                const double b_sym = 0.5 * (bij + bji);

                const double u_pair = fcij * (VR - b_sym * VA);
                u[i] += 0.5 * u_pair;
                u[j] += 0.5 * u_pair;

                const double denom_ij = one_p_Xij * sqrt_one_p_Xij;
                const double denom_ji = one_p_Xji * sqrt_one_p_Xji;
                const double fk_ij = (denom_ij > 0.0) ? (fcij * VA * (0.25 / denom_ij)) : 0.0;
                const double fk_ji = (denom_ji > 0.0) ? (fcij * VA * (0.25 / denom_ji)) : 0.0;

                for (int idx = 0; idx < count_ij; ++idx)
                {
                    const int kk = buf_ij.k_list[idx];
                    fx[i] += -fk_ij * buf_ij.d_center_x[idx];
                    fy[i] += -fk_ij * buf_ij.d_center_y[idx];
                    fz[i] += -fk_ij * buf_ij.d_center_z[idx];

                    fx[j] += -fk_ij * buf_ij.d_other_x[idx];
                    fy[j] += -fk_ij * buf_ij.d_other_y[idx];
                    fz[j] += -fk_ij * buf_ij.d_other_z[idx];

                    fx[kk] += -fk_ij * buf_ij.d_k_x[idx];
                    fy[kk] += -fk_ij * buf_ij.d_k_y[idx];
                    fz[kk] += -fk_ij * buf_ij.d_k_z[idx];
                }

                for (int idx = 0; idx < count_ji; ++idx)
                {
                    const int kk = buf_ji.k_list[idx];
                    fx[j] += -fk_ji * buf_ji.d_center_x[idx];
                    fy[j] += -fk_ji * buf_ji.d_center_y[idx];
                    fz[j] += -fk_ji * buf_ji.d_center_z[idx];

                    fx[i] += -fk_ji * buf_ji.d_other_x[idx];
                    fy[i] += -fk_ji * buf_ji.d_other_y[idx];
                    fz[i] += -fk_ji * buf_ji.d_other_z[idx];

                    fx[kk] += -fk_ji * buf_ji.d_k_x[idx];
                    fy[kk] += -fk_ji * buf_ji.d_k_y[idx];
                    fz[kk] += -fk_ji * buf_ji.d_k_z[idx];
                }

                const double Fij_scalar = nfcij * (VR - b_sym * VA) + fcij * (nVR - b_sym * nVA);

                fx[i] += seij_x * Fij_scalar;
                fy[i] += seij_y * Fij_scalar;
                fz[i] += seij_z * Fij_scalar;
                fx[j] -= seij_x * Fij_scalar;
                fy[j] -= seij_y * Fij_scalar;
                fz[j] -= seij_z * Fij_scalar;
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        data.atoms[i].f.a00 = fx[i];
        data.atoms[i].f.a10 = fy[i];
        data.atoms[i].f.a20 = fz[i];
        U_atom[i] = u[i];

        data.F_all.a00 += fx[i];
        data.F_all.a10 += fy[i];
        data.F_all.a20 += fz[i];
    }
}
