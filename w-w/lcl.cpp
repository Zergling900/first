#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <omp.h>

#include "3.h"
#include "void.h"

void lcl0(Data &data, Cell_List &cl, const parameter1 &pr1, const parameter2 &pr2)
{
    int n = data.n;
    double cutoff = pr2.R + pr2.D;
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
        std::vector<int> k_list_ij;
        std::vector<int> k_list_ji;

        std::vector<Matrix31> dXij_dri_list;
        std::vector<Matrix31> dXij_drj_list;
        std::vector<Matrix31> dXij_drk_list;

        std::vector<Matrix31> dXji_dri_list;
        std::vector<Matrix31> dXji_drj_list;
        std::vector<Matrix31> dXji_drk_list;

        k_list_ij.reserve(n);
        k_list_ji.reserve(n);
        dXij_dri_list.reserve(n);
        dXij_drj_list.reserve(n);
        dXij_drk_list.reserve(n);
        dXji_dri_list.reserve(n);
        dXji_drj_list.reserve(n);
        dXji_drk_list.reserve(n);

#pragma omp for schedule(guided, 32) reduction(+:fx[:n], fy[:n], fz[:n], u[:n])
        for (int i = 0; i < n; i++)
        {
            const Atom &ai = data.atoms[i];
            const std::string &ni = ai.name;

            int c_i = cl.Cell[i];
            int cz_i = c_i / (Mx * My);
            int cc_i = c_i - cz_i * (Mx * My); // just a number
            int cy_i = cc_i / Mx;
            int cx_i = cc_i - cy_i * Mx;
            // #pragma omp parallel for schedule(guided, 32)
            for (int dcz = -1; dcz <= 1; dcz++)
            {
                for (int dcy = -1; dcy <= 1; dcy++)
                {
                    for (int dcx = -1; dcx <= 1; dcx++)
                    {
                        int cz2 = wrap_cell(cz_i + dcz, Mz);
                        int cy2 = wrap_cell(cy_i + dcy, My);
                        int cx2 = wrap_cell(cx_i + dcx, Mx);
                        int nc = cx2 + cy2 * Mx + cz2 * Mx * My;

                        int begin = cl.cell_offset[nc];
                        int end = cl.cell_offset[nc + 1];

                        for (int jj = begin; jj < end; jj++)
                        {
                            int j = cl.atom_indices[jj];
                            if (j <= i)
                            {
                                continue;
                            }

                            const Atom &aj = data.atoms[j];
                            const std::string &nj = aj.name;

                            // Choose ABOP parameter set based on species of i and j
                            const parameter2 *p = nullptr;

                            bool i_is_W = (ni == "W");
                            bool j_is_W = (nj == "W");
                            bool i_is_Be = (ni == "Be");
                            bool j_is_Be = (nj == "Be");

                            if (i_is_W && j_is_W)
                            {
                                p = &pr2_WW;
                            }
                            else if (i_is_Be && j_is_Be)
                            {
                                p = &pr2_BB;
                            }
                            else
                            {
                                p = &pr2_WB;
                            }

                            double D0 = p->D0;
                            double R = p->R;
                            double D = p->D;
                            double mu = p->mu; // 2mu in paper
                            double r0 = p->r0;
                            double beta = p->beta;
                            double S = p->S;
                            double gamma = p->gamma;

                            // !!! rename to avoid shadowing cell index c_i / etc.
                            double c_ang = p->c;
                            double d_ang = p->d;
                            double h_ang = p->h;

                            // displacement r_ij = r_j - r_i
                            Matrix31 drij = aj.r - ai.r;

                            // PBC (real space)
                            if (drij.a00 < -Lxh)
                                drij.a00 += Lx;
                            else if (drij.a00 >= Lxh)
                                drij.a00 -= Lx;

                            if (drij.a10 < -Lyh)
                                drij.a10 += Ly;
                            else if (drij.a10 >= Lyh)
                                drij.a10 -= Ly;

                            if (drij.a20 < -Lzh)
                                drij.a20 += Lz;
                            else if (drij.a20 >= Lzh)
                                drij.a20 -= Lz;

                            double rij2 = drij.a00 * drij.a00 + drij.a10 * drij.a10 + drij.a20 * drij.a20;
                            if (rij2 == 0.0)
                            {
                                continue;
                            }

                            double rij = std::sqrt(rij2);
                            Matrix31 seij = drij * (1.0 / rij);
                            Matrix31 drji = seij * (-rij);
                            double rji = rij;
                            Matrix31 seji = seij * (-1.0);

                            // Pair VR / VA
                            double bsR = beta * std::sqrt(2.0 * S);
                            double bsA = beta * std::sqrt(2.0 / S);

                            double VR = (D0 / (S - 1.0)) * std::exp(-bsR * (rij - r0));
                            double VA = (S * D0 / (S - 1.0)) * std::exp(-bsA * (rij - r0));

                            // Cutoff fc_ij and derivative
                            double fcij, nfcij;


                            double inside = (rij > (R - D)) && (rij < (R + D));
                            double below = (rij <= (R - D));
                            double above = (rij >= (R + D));
                            double x = 0.5 * M_PI * (rij - R) / D;
                            fcij = below * 1.0 + inside * (0.5 - 0.5 * sin(x)) + above * 0.0;
                            nfcij = inside * (-0.5 * cos(x) * (0.5 * M_PI / D));

                            if (fcij == 0.0)
                            {
                                continue;
                            }

                            double nVR = (D0 / (S - 1.0)) * (-bsR) * std::exp(-bsR * (rij - r0));
                            double nVA = (S * D0 / (S - 1.0)) * (-bsA) * std::exp(-bsA * (rij - r0));

                            // -------------------------------------------------
                            // Three-body: build Xij using k in i-neighborhood,
                            //            build Xji using k in j-neighborhood
                            // -------------------------------------------------
                            double Xij = 0.0;
                            double Xji = 0.0;

                            k_list_ij.clear();
                            k_list_ji.clear();

                            dXij_dri_list.clear();
                            dXij_drj_list.clear();
                            dXij_drk_list.clear();

                            dXji_dri_list.clear();
                            dXji_drj_list.clear();
                            dXji_drk_list.clear();
                            auto accumulate_X_center = [&](int center,
                                                           const Atom &ac,
                                                           const Matrix31 &dr_co,
                                                           const Matrix31 &se_co,
                                                           double r_co,
                                                           double &X,
                                                           std::vector<Matrix31> &dX_dr_center_list,
                                                           std::vector<Matrix31> &dX_dr_other_list,
                                                           std::vector<Matrix31> &dX_drk_list,
                                                           std::vector<int> &k_list)
                            {
                                for_each_in_27_cells(center, [&](int k)
                                                     {//function:
                                    if (k == i || k == j)
                                    {
                                        return;
                                    }

                                    const Atom &ak = data.atoms[k];
                                    Matrix31 dr_ck = ak.r - ac.r;

                                    if (dr_ck.a00 < -Lxh)        dr_ck.a00 += Lx;
                                    else if (dr_ck.a00 >= Lxh)   dr_ck.a00 -= Lx;

                                    if (dr_ck.a10 < -Lyh)        dr_ck.a10 += Ly;
                                    else if (dr_ck.a10 >= Lyh)   dr_ck.a10 -= Ly;

                                    if (dr_ck.a20 < -Lzh)        dr_ck.a20 += Lz;
                                    else if (dr_ck.a20 >= Lzh)   dr_ck.a20 -= Lz;

                                    double r_ck2 = dr_ck.a00 * dr_ck.a00
                                                 + dr_ck.a10 * dr_ck.a10
                                                 + dr_ck.a20 * dr_ck.a20;
                                    if (r_ck2 == 0.0)
                                    {
                                        return;
                                    }
                                    double r_ck = std::sqrt(r_ck2);

                                    double fc_ck, dfc_dr_ck;
                                    if (r_ck <= R - D)
                                    {
                                        fc_ck = 1.0;
                                        dfc_dr_ck = 0.0;
                                    }
                                    else if (r_ck >= R + D)
                                    {
                                        fc_ck = 0.0;
                                        dfc_dr_ck = 0.0;
                                    }
                                    else
                                    {
                                        double xx = 0.5 * M_PI * (r_ck - R) / D;
                                        fc_ck = 0.5 - 0.5 * std::sin(xx);
                                        dfc_dr_ck = -0.5 * std::cos(xx) * (0.5 * M_PI / D);
                                    }

                                    if (fc_ck <= 0.0)
                                    {
                                        return;
                                    }

                                    Matrix31 se_ck = dr_ck * (1.0 / r_ck);

                                    double dot_co_ck = dr_co.a00 * dr_ck.a00 + dr_co.a10 * dr_ck.a10 + dr_co.a20 * dr_ck.a20;
                                    double cost_cok = dot_co_ck / (r_co * r_ck);
                                    if (cost_cok > 1.0) cost_cok = 1.0;
                                    if (cost_cok < -1.0) cost_cok = -1.0;

                                    double u_ang = d_ang * d_ang + (h_ang + cost_cok) * (h_ang + cost_cok);
                                    double g_ck = gamma * (1.0 + c_ang * c_ang * (1.0 / (d_ang * d_ang) - 1.0 / u_ang));
                                    double dg_dcos = gamma * c_ang * c_ang * 2.0 * (h_ang + cost_cok) / (u_ang * u_ang);

                                    double inv_r_co = 1.0 / r_co;
                                    double inv_r_ck = 1.0 / r_ck;
                                    double inv_r_co2 = inv_r_co * inv_r_co;
                                    double inv_r_ck2 = inv_r_ck * inv_r_ck;
                                    double inv_r_co_ck = inv_r_co * inv_r_ck;

                                    Matrix31 term_dr_co = dr_co * (dot_co_ck * inv_r_co2);
                                    Matrix31 dcos_ddr_co = (dr_ck - term_dr_co) * inv_r_co_ck;

                                    Matrix31 term_dr_ck = dr_ck * (dot_co_ck * inv_r_ck2);
                                    Matrix31 dcos_ddr_ck = (dr_co - term_dr_ck) * inv_r_co_ck;

                                    Matrix31 dcos_drc = (dcos_ddr_co + dcos_ddr_ck) * (-1.0);
                                    Matrix31 dcos_dro = dcos_ddr_co;
                                    Matrix31 dcos_drk = dcos_ddr_ck;

                                    Matrix31 dg_drc = dcos_drc * dg_dcos;
                                    Matrix31 dg_dro = dcos_dro * dg_dcos;
                                    Matrix31 dg_drk = dcos_drk * dg_dcos;

                                    double e_cok = std::exp(mu * (r_co - r_ck));

                                    Matrix31 dfc_drc = se_ck * (-dfc_dr_ck);
                                    Matrix31 dfc_dro(0.0, 0.0, 0.0);
                                    Matrix31 dfc_drk = se_ck * dfc_dr_ck;

                                    Matrix31 de_drc = (se_ck - se_co) * (mu * e_cok);
                                    Matrix31 de_dro = se_co * (mu * e_cok);
                                    Matrix31 de_drk = se_ck * (-mu * e_cok);

                                    Matrix31 dX_drc = dfc_drc * (g_ck * e_cok) + dg_drc * (fc_ck * e_cok) + de_drc * (fc_ck * g_ck);
                                    Matrix31 dX_dro = dfc_dro * (g_ck * e_cok) + dg_dro * (fc_ck * e_cok) + de_dro * (fc_ck * g_ck);
                                    Matrix31 dX_drk = dfc_drk * (g_ck * e_cok) + dg_drk * (fc_ck * e_cok) + de_drk * (fc_ck * g_ck);
                                    X += fc_ck * g_ck * e_cok;

                                    dX_dr_center_list.push_back(dX_drc);
                                    dX_dr_other_list.push_back(dX_dro);
                                    dX_drk_list.push_back(dX_drk);
                                    k_list.push_back(k);
                                    //function_end
                                });
                            };

                            // --------------------
                            // Xij : i-centered
                            // --------------------
                            accumulate_X_center(i, ai, drij, seij, rij, Xij,
                                                dXij_dri_list, dXij_drj_list, dXij_drk_list, k_list_ij);

                            // --------------------
                            // Xji : j-centered
                            // --------------------
                            accumulate_X_center(j, aj, drji, seji, rji, Xji,
                                                dXji_drj_list, dXji_dri_list, dXji_drk_list, k_list_ji);

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

                                Matrix31 dXij_dri = dXij_dri_list[idx];
                                Matrix31 dXij_drj = dXij_drj_list[idx];
                                Matrix31 dXij_drk = dXij_drk_list[idx];

                                Matrix31 Fi_xyz = dXij_dri * (-fk_ij);
                                Matrix31 Fj_xyz = dXij_drj * (-fk_ij);
                                Matrix31 Fk_xyz = dXij_drk * (-fk_ij);

                                fx[i] += Fi_xyz.a00;
                                fy[i] += Fi_xyz.a10;
                                fz[i] += Fi_xyz.a20;

                                fx[j] += Fj_xyz.a00;
                                fy[j] += Fj_xyz.a10;
                                fz[j] += Fj_xyz.a20;

                                fx[kk] += Fk_xyz.a00;
                                fy[kk] += Fk_xyz.a10;
                                fz[kk] += Fk_xyz.a20;
                            }

                            // forces from Xji (j-centered)
                            for (size_t idx = 0; idx < k_list_ji.size(); ++idx)
                            {
                                int kk = k_list_ji[idx];

                                Matrix31 dXji_dri = dXji_dri_list[idx];
                                Matrix31 dXji_drj = dXji_drj_list[idx];
                                Matrix31 dXji_drk = dXji_drk_list[idx];

                                Matrix31 Fi_xyz = dXji_dri * (-fk_ji);
                                Matrix31 Fj_xyz = dXji_drj * (-fk_ji);
                                Matrix31 Fk_xyz = dXji_drk * (-fk_ji);

                                fx[i] += Fi_xyz.a00;
                                fy[i] += Fi_xyz.a10;
                                fz[i] += Fi_xyz.a20;

                                fx[j] += Fj_xyz.a00;
                                fy[j] += Fj_xyz.a10;
                                fz[j] += Fj_xyz.a20;

                                fx[kk] += Fk_xyz.a00;
                                fy[kk] += Fk_xyz.a10;
                                fz[kk] += Fk_xyz.a20;
                            }
                            // pair (radial) forces from fc, VR, VA
                            double Fij_scalar = nfcij * (VR - b_sym * VA) + fcij * (nVR - b_sym * nVA);

                            Matrix31 fij_xyz = seij * Fij_scalar;
                            fx[i] += fij_xyz.a00;
                            fy[i] += fij_xyz.a10;
                            fz[i] += fij_xyz.a20;
                            fx[j] -= fij_xyz.a00;
                            fy[j] -= fij_xyz.a10;
                            fz[j] -= fij_xyz.a20;
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
