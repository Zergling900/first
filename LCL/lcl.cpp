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

    if (cl.Mx < 1) cl.Mx = 1;
    if (cl.My < 1) cl.My = 1;
    if (cl.Mz < 1) cl.Mz = 1;

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
    #pragma omp parallel for schedule(guided, 32)
    for (int i = 0; i < n; i++)
    {
        const Atom &ai = data.atoms[i];
        const std::string &ni = ai.name;

        int c_i = cl.Cell[i];
        int cz_i = c_i / (Mx * My);
        int cc_i = c_i - cz_i * (Mx * My); // just a number
        int cy_i = cc_i / Mx;
        int cx_i = cc_i - cy_i * Mx;

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
/*
                        if (rij <= (R - D))
                        {
                            fcij = 1.0;
                            nfcij = 0.0;
                        }
                        else if (rij >= (R + D))
                        {
                            fcij = 0.0;
                            nfcij = 0.0;
                        }
                        else
                        {
                            double x = 0.5 * M_PI * (rij - R) / D;
                            fcij = 0.5 - 0.5 * std::sin(x);
                            nfcij = -0.5 * std::cos(x) * (0.5 * M_PI / D);
                        }
*/

                        double inside = (rij > (R - D)) && (rij < (R + D));
                        double below = (rij <= (R - D)); double above = (rij >= (R + D));
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

                        // --------------------
                        // Xij : i-centered
                        // --------------------
                        for_each_in_27_cells(i, [&](int k)
                                             {
                            if (k == i || k == j)
                            {
                                return;
                            }

                            const Atom &ak = data.atoms[k];

                            Matrix31 drik = ak.r - ai.r;

                            if (drik.a00 < -Lxh)        drik.a00 += Lx;
                            else if (drik.a00 >= Lxh)   drik.a00 -= Lx;

                            if (drik.a10 < -Lyh)        drik.a10 += Ly;
                            else if (drik.a10 >= Lyh)   drik.a10 -= Ly;

                            if (drik.a20 < -Lzh)        drik.a20 += Lz;
                            else if (drik.a20 >= Lzh)   drik.a20 -= Lz;

                            double rik2 = drik.a00 * drik.a00
                                        + drik.a10 * drik.a10
                                        + drik.a20 * drik.a20;
                            if (rik2 == 0.0)
                            {
                                return;
                            }
                            double rik = std::sqrt(rik2);

                            double fcik, dfc_drik;
                            if (rik <= R - D)
                            {
                                fcik = 1.0;
                                dfc_drik = 0.0;
                            }
                            else if (rik >= R + D)
                            {
                                fcik = 0.0;
                                dfc_drik = 0.0;
                            }
                            else
                            {
                                double xx = 0.5 * M_PI * (rik - R) / D;
                                fcik = 0.5 - 0.5 * std::sin(xx);
                                dfc_drik = -0.5 * std::cos(xx) * (0.5 * M_PI / D);
                            }

                            if (fcik <= 0.0)
                            {
                                return;
                            }

                            Matrix31 dXij_dri(0.0, 0.0, 0.0);
                            Matrix31 dXij_drj(0.0, 0.0, 0.0);
                            Matrix31 dXij_drk(0.0, 0.0, 0.0);

                            Matrix31 seik = drik * (1.0 / rik);

                            double dot_ij_ik = drij.a00 * drik.a00 + drij.a10 * drik.a10 + drij.a20 * drik.a20;
                            double cost_ijk = dot_ij_ik / (rij * rik);
                            if (cost_ijk > 1.0) cost_ijk = 1.0;
                            if (cost_ijk < -1.0) cost_ijk = -1.0;

                            double u_ang = d_ang * d_ang + (h_ang + cost_ijk) * (h_ang + cost_ijk);
                            double gik = gamma * (1.0 + c_ang * c_ang * (1.0 / (d_ang * d_ang) - 1.0 / u_ang));
                            double dg_dcos = gamma * c_ang * c_ang * 2.0 * (h_ang + cost_ijk) / (u_ang * u_ang);

                            double inv_rij = 1.0 / rij;
                            double inv_rik = 1.0 / rik;
                            double inv_rij2 = inv_rij * inv_rij;
                            double inv_rik2 = inv_rik * inv_rik;
                            double inv_rij_rik = inv_rij * inv_rik;

                            Matrix31 term_drij = drij * (dot_ij_ik * inv_rij2);
                            Matrix31 dcos_ddrij = (drik - term_drij) * inv_rij_rik;

                            Matrix31 term_drik = drik * (dot_ij_ik * inv_rik2);
                            Matrix31 dcos_ddrik = (drij - term_drik) * inv_rij_rik;

                            Matrix31 dcos_dri = (dcos_ddrij + dcos_ddrik) * (-1.0);
                            Matrix31 dcos_drj = dcos_ddrij;
                            Matrix31 dcos_drk = dcos_ddrik;

                            Matrix31 dg_dri = dcos_dri * dg_dcos;
                            Matrix31 dg_drj = dcos_drj * dg_dcos;
                            Matrix31 dg_drk = dcos_drk * dg_dcos;

                            double eijk = std::exp(mu * (rij - rik));

                            Matrix31 dfc_dri = seik * (-dfc_drik);
                            Matrix31 dfc_drj(0.0, 0.0, 0.0);
                            Matrix31 dfc_drk = seik * dfc_drik;

                            Matrix31 de_dri = (seik - seij) * (mu * eijk);
                            Matrix31 de_drj = seij * (mu * eijk);
                            Matrix31 de_drk = seik * (-mu * eijk);

                            Matrix31 dX_dri = dfc_dri * (gik * eijk) + dg_dri * (fcik * eijk) + de_dri * (fcik * gik);
                            Matrix31 dX_drj = dfc_drj * (gik * eijk) + dg_drj * (fcik * eijk) + de_drj * (fcik * gik);
                            Matrix31 dX_drk = dfc_drk * (gik * eijk) + dg_drk * (fcik * eijk) + de_drk * (fcik * gik);

                            Xij += fcik * gik * eijk;

                            dXij_dri = dX_dri;
                            dXij_drj = dX_drj;
                            dXij_drk = dX_drk;

                            dXij_dri_list.push_back(dXij_dri);
                            dXij_drj_list.push_back(dXij_drj);
                            dXij_drk_list.push_back(dXij_drk);
                            k_list_ij.push_back(k); });

                        // --------------------
                        // Xji : j-centered
                        // --------------------
                        for_each_in_27_cells(j, [&](int k)
                                             {
                            if (k == i || k == j)
                            {
                                return;
                            }

                            const Atom &ak = data.atoms[k];

                            Matrix31 drjk = ak.r - aj.r;

                            if (drjk.a00 < -Lxh) drjk.a00 += Lx;
                            else if (drjk.a00 >= Lxh) drjk.a00 -= Lx;

                            if (drjk.a10 < -Lyh) drjk.a10 += Ly;
                            else if (drjk.a10 >= Lyh) drjk.a10 -= Ly;

                            if (drjk.a20 < -Lzh) drjk.a20 += Lz;
                            else if (drjk.a20 >= Lzh) drjk.a20 -= Lz;

                            double rjk2 = drjk.a00 * drjk.a00 + drjk.a10 * drjk.a10 + drjk.a20 * drjk.a20;
                            if (rjk2 == 0.0)
                            {
                                return;
                            }
                            double rjk = std::sqrt(rjk2);

                            double fcjk, dfc_drjk;
                            if (rjk <= R - D)
                            {
                                fcjk = 1.0;
                                dfc_drjk = 0.0;
                            }
                            else if (rjk >= R + D)
                            {
                                fcjk = 0.0;
                                dfc_drjk = 0.0;
                            }
                            else
                            {
                                double xx = 0.5 * M_PI * (rjk - R) / D;
                                fcjk = 0.5 - 0.5 * std::sin(xx);
                                dfc_drjk = -0.5 * std::cos(xx) * (0.5 * M_PI / D);
                            }

                            if (fcjk <= 0.0)
                            {
                                return;
                            }

                            Matrix31 dXji_dri(0.0, 0.0, 0.0);
                            Matrix31 dXji_drj(0.0, 0.0, 0.0);
                            Matrix31 dXji_drk(0.0, 0.0, 0.0);

                            Matrix31 sejk = drjk * (1.0 / rjk);

                            double dot_ji_jk = drji.a00 * drjk.a00 + drji.a10 * drjk.a10 + drji.a20 * drjk.a20;
                            double cost_jik = dot_ji_jk / (rji * rjk);
                            if (cost_jik > 1.0) cost_jik = 1.0;
                            if (cost_jik < -1.0) cost_jik = -1.0;

                            double u_ang2 = d_ang * d_ang + (h_ang + cost_jik) * (h_ang + cost_jik);
                            double gjk = gamma * (1.0 + c_ang * c_ang * (1.0 / (d_ang * d_ang) - 1.0 / u_ang2));
                            double dg2_dcos = gamma * c_ang * c_ang * 2.0 * (h_ang + cost_jik) / (u_ang2 * u_ang2);

                            double inv_rji = 1.0 / rji;
                            double inv_rjk = 1.0 / rjk;
                            double inv_rji2 = inv_rji * inv_rji;
                            double inv_rjk2 = inv_rjk * inv_rjk;
                            double inv_rji_rjk = inv_rji * inv_rjk;

                            Matrix31 term_drji = drji * (dot_ji_jk * inv_rji2);
                            Matrix31 dcos_ddrji = (drjk - term_drji) * inv_rji_rjk;

                            Matrix31 term_drjk = drjk * (dot_ji_jk * inv_rjk2);
                            Matrix31 dcos_ddrjk = (drji - term_drjk) * inv_rji_rjk;

                            Matrix31 dcos_drj = (dcos_ddrji + dcos_ddrjk) * (-1.0);
                            Matrix31 dcos_dri = dcos_ddrji;
                            Matrix31 dcos_drk = dcos_ddrjk;

                            Matrix31 dg_dri2 = dcos_dri * dg2_dcos;
                            Matrix31 dg_drj2 = dcos_drj * dg2_dcos;
                            Matrix31 dg_drk2 = dcos_drk * dg2_dcos;

                            double ejik = std::exp(mu * (rji - rjk));

                            Matrix31 dfc_dri2(0.0, 0.0, 0.0);
                            Matrix31 dfc_drj2 = sejk * (-dfc_drjk);
                            Matrix31 dfc_drk2 = sejk * dfc_drjk;

                            Matrix31 de_drj2 = (sejk - seji) * (mu * ejik);
                            Matrix31 de_dri2 = seji * (mu * ejik);
                            Matrix31 de_drk2 = sejk * (-mu * ejik);

                            Matrix31 dX2_dri = dfc_dri2 * (gjk * ejik) + dg_dri2 * (fcjk * ejik) + de_dri2 * (fcjk * gjk);
                            Matrix31 dX2_drj = dfc_drj2 * (gjk * ejik) + dg_drj2 * (fcjk * ejik) + de_drj2 * (fcjk * gjk);
                            Matrix31 dX2_drk = dfc_drk2 * (gjk * ejik) + dg_drk2 * (fcjk * ejik) + de_drk2 * (fcjk * gjk);

                            Xji += fcjk * gjk * ejik;

                            dXji_dri = dX2_dri;
                            dXji_drj = dX2_drj;
                            dXji_drk = dX2_drk;

                            dXji_dri_list.push_back(dXji_dri);
                            dXji_drj_list.push_back(dXji_drj);
                            dXji_drk_list.push_back(dXji_drk);
                            k_list_ji.push_back(k); });

                        // bond order terms
                        double bij = 1.0 / std::sqrt(1.0 + Xij);
                        double bji = 1.0 / std::sqrt(1.0 + Xji);
                        double b_sym = 0.5 * (bij + bji);

                        double u_pair = fcij * (VR - b_sym * VA);

                        U_atom[i] += 0.5 * u_pair;
                        U_atom[j] += 0.5 * u_pair;

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

                            data.atoms[i].f = data.atoms[i].f + Fi_xyz;
                            data.atoms[j].f = data.atoms[j].f + Fj_xyz;
                            data.atoms[kk].f = data.atoms[kk].f + Fk_xyz;
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

                            data.atoms[i].f = data.atoms[i].f + Fi_xyz;
                            data.atoms[j].f = data.atoms[j].f + Fj_xyz;
                            data.atoms[kk].f = data.atoms[kk].f + Fk_xyz;
                        }

                        // pair (radial) forces from fc, VR, VA
                        double Fij_scalar = nfcij * (VR - b_sym * VA) + fcij * (nVR - b_sym * nVA);

                        Matrix31 fij_xyz = seij * Fij_scalar;

                        data.atoms[i].f = data.atoms[i].f + fij_xyz;
                        data.atoms[j].f = data.atoms[j].f - fij_xyz;
                    }
                }
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        data.F_all = data.F_all + data.atoms[i].f;
    }
}

void lcl222(Data &data, Cell_List &cl, const parameter1 &pr1, const parameter2 &pr2, vector<double> &U_atom)
{
    int n = data.n;
    double Wx = cl.Wx;
    double Wy = cl.Wy;
    double Wz = cl.Wz;

    int Mx = cl.Mx;
    int My = cl.My;
    int Mz = cl.Mz;
    for (int i = 0; i < n; i++)
    {
        int c = cl.Cell[i];
        int cz = c / (Mx * My);
        int cc = c - cz * (Mx * My); // just a number
        int cy = cc / Mx;
        int cx = cc - cy * Mx;

        int cz2, cy2, cx2;
        for (int dcz = -1; dcz <= 1; dcz++)
        {
            for (int dcy = -1; dcy <= 1; dcy++)
            {
                for (int dcx = -1; dcx <= 1; dcx++)
                {
                    int cz2 = cz + dcz;
                    int cy2 = cy + dcy;
                    int cx2 = cx + dcx;

                    // CELL PBC
                    if (cx2 < 0)
                    {
                        cx2 += Mx;
                    }
                    if (cx2 >= Mx)
                    {
                        cx2 -= Mx;
                    }
                    if (cy2 < 0)
                    {
                        cy2 += My;
                    }
                    if (cy2 >= My)
                    {
                        cy2 -= My;
                    }
                    if (cz2 < 0)
                    {
                        cz2 += Mz;
                    }
                    if (cz2 >= Mz)
                    {
                        cz2 -= Mz;
                    }

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
                    }
                }
            }
        }
    }
}
