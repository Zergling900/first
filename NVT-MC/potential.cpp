#include <iostream>
// #include <random>
#include <sstream>
#include <string>
#include <cmath>
// #include <fstream>
#include <vector>
// #include <iomanip>

#include "3.h"

using namespace std;

// ----------------------------
// potential
// ----------------------------

// LJ potential U = 4*epsilon*(sigma^12/r^12 - sigma^6/r^6)

void BeW_potential2(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,
                   vector<double> &U_atom)
{
    int i, j, k;

    // Reset total force and per-atom forces
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    U_atom.assign(data.n, 0.0);

    double Lx  = data.Box.a00;
    double Ly  = data.Box.a11;
    double Lz  = data.Box.a22;
    double Lxh = Lx / 2.0;
    double Lyh = Ly / 2.0;
    double Lzh = Lz / 2.0;

    for (i = 0; i < data.n; i++)
    {
        data.atoms[i].f = Matrix31(0.0, 0.0, 0.0);
    }

    // i < j double loop so that each unordered pair is handled exactly once
    for (i = 0; i < data.n; i++)
    {
        const Atom &ai = data.atoms[i];
        const std::string &ni = ai.name;

        for (j = i + 1; j < data.n; j++)
        {
            const Atom &aj = data.atoms[j];
            const std::string &nj = aj.name;

            // Choose ABOP parameter set based on species of i and j
            const parameter2 *p = nullptr;

            bool i_is_W  = (ni == "W");
            bool j_is_W  = (nj == "W");
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
                // mixed W-Be or Be-W pair
                p = &pr2_WB;
            }

            // Extract ABOP parameters for this i-j pair
            double D0    = p->D0;
            double R     = p->R;
            double D     = p->D;
            double mu    = p->mu;   // here mu corresponds to "2mu" in the paper
            double r0    = p->r0;
            double beta  = p->beta;
            double S     = p->S;
            double gamma = p->gamma;
            double c     = p->c;
            double d     = p->d;
            double h     = p->h;
            // double rf  = p->rf; // !ZBL-related, not used yet
            // double bf  = p->bf;

            // displacement r_ij = r_j - r_i
            Matrix31 drij = aj.r - ai.r;
            //PBC*******************************************
            if (drij.a00 < -Lxh)      drij.a00 += Lx;
            else if (drij.a00 >= Lxh) drij.a00 -= Lx;

            if (drij.a10 < -Lyh)      drij.a10 += Ly;
            else if (drij.a10 >= Lyh) drij.a10 -= Ly;

            if (drij.a20 < -Lzh)      drij.a20 += Lz;
            else if (drij.a20 >= Lzh) drij.a20 -= Lz;
            //PBC*******************************************
            double rij2 = drij.a00 * drij.a00
                        + drij.a10 * drij.a10
                        + drij.a20 * drij.a20;
            if (rij2 == 0.0)
                continue;

            double rij = std::sqrt(rij2);
            Matrix31 seij = drij * (1.0 / rij);    // unit vector from i to j
            Matrix31 drji = seij * (-rij);         // r_ji = -r_ij, same magnitude
            double  rji   = rij;
            Matrix31 seji = seij * (-1.0);

            // Pair VR / VA
            double bsR = beta * std::sqrt(2.0 * S);
            double bsA = beta * std::sqrt(2.0 / S);

            double VR  = (D0 / (S - 1.0)) * std::exp(-bsR * (rij - r0));
            double VA  = (S * D0 / (S - 1.0)) * std::exp(-bsA * (rij - r0));

            // Cutoff fc_ij and its radial derivative nfcij = dfc/dr
            double fcij, nfcij;

            if (rij <= (R - D))
            {
                fcij  = 1.0;
                nfcij = 0.0;
            }
            else if (rij >= (R + D))
            {
                fcij  = 0.0;
                nfcij = 0.0;
            }
            else
            {
                double x = 0.5 * M_PI * (rij - R) / D;
                fcij  = 0.5 - 0.5 * std::sin(x);
                nfcij = -0.5 * std::cos(x) * (0.5 * M_PI / D);
            }

            // If this pair is fully outside the cutoff, skip all k-terms
            if (fcij == 0.0)
            {
                continue;
            }

            
            double nVR = (D0 / (S - 1.0)) * (-bsR) * std::exp(-bsR * (rij - r0));
            double nVA = (S * D0 / (S - 1.0)) * (-bsA) * std::exp(-bsA * (rij - r0));

            // Three-body: Xij and Xji plus their gradients
            double Xij = 0.0;
            double Xji = 0.0;

            std::vector<Matrix31> dXij_dri_list;
            std::vector<Matrix31> dXij_drj_list;
            std::vector<Matrix31> dXij_drk_list;

            std::vector<Matrix31> dXji_dri_list;
            std::vector<Matrix31> dXji_drj_list;
            std::vector<Matrix31> dXji_drk_list;

            std::vector<int> k_list;

            dXij_dri_list.reserve(data.n);
            dXij_drj_list.reserve(data.n);
            dXij_drk_list.reserve(data.n);
            dXji_dri_list.reserve(data.n);
            dXji_drj_list.reserve(data.n);
            dXji_drk_list.reserve(data.n);
            k_list.reserve(data.n);

            for (k = 0; k < data.n; k++)
            {
                if (k == i || k == j)
                    continue;

                const Atom &ak = data.atoms[k];

                // ---- central i: vectors r_ik ----
                Matrix31 drik = ak.r - ai.r;
                // PBC*******************************************
                if (drik.a00 < -Lxh)      drik.a00 += Lx;
                else if (drik.a00 >= Lxh) drik.a00 -= Lx;

                if (drik.a10 < -Lyh)      drik.a10 += Ly;
                else if (drik.a10 >= Lyh) drik.a10 -= Ly;

                if (drik.a20 < -Lzh)      drik.a20 += Lz;
                else if (drik.a20 >= Lzh) drik.a20 -= Lz;
                // PBC*******************************************

                double rik2 = drik.a00 * drik.a00
                            + drik.a10 * drik.a10
                            + drik.a20 * drik.a20;
                double rik  = std::sqrt(rik2);

                // ---- central j: vectors r_jk ----
                Matrix31 drjk = ak.r - aj.r;
                if (drjk.a00 < -Lxh)      drjk.a00 += Lx;
                else if (drjk.a00 >= Lxh) drjk.a00 -= Lx;

                if (drjk.a10 < -Lyh)      drjk.a10 += Ly;
                else if (drjk.a10 >= Lyh) drjk.a10 -= Ly;

                if (drjk.a20 < -Lzh)      drjk.a20 += Lz;
                else if (drjk.a20 >= Lzh) drjk.a20 -= Lz;

                double rjk2 = drjk.a00 * drjk.a00
                            + drjk.a10 * drjk.a10
                            + drjk.a20 * drjk.a20;
                double rjk  = std::sqrt(rjk2);

                // For this k, initialize gradient accumulators (for Xij and Xji)
                Matrix31 dXij_dri(0.0, 0.0, 0.0);
                Matrix31 dXij_drj(0.0, 0.0, 0.0);
                Matrix31 dXij_drk(0.0, 0.0, 0.0);

                Matrix31 dXji_dri(0.0, 0.0, 0.0);
                Matrix31 dXji_drj(0.0, 0.0, 0.0);
                Matrix31 dXji_drk(0.0, 0.0, 0.0);

                // ---------- Xij part ----------
                if (rik > 0.0)
                {
                    // cutoff fc(r_ik)
                    double fcik, dfc_drik;
                    if (rik <= R - D)
                    {
                        fcik     = 1.0;
                        dfc_drik = 0.0;
                    }
                    else if (rik >= R + D)
                    {
                        fcik     = 0.0;
                        dfc_drik = 0.0;
                    }
                    else
                    {
                        double xx = 0.5 * M_PI * (rik - R) / D;
                        fcik     = 0.5 - 0.5 * std::sin(xx);
                        dfc_drik = -0.5 * std::cos(xx) * (0.5 * M_PI / D);
                    }

                    if (fcik > 0.0)
                    {
                        Matrix31 seik = drik * (1.0 / rik);

                        // cos(theta_ijk) with central i: angle between r_ij and r_ik
                        double dot_ij_ik = drij.a00 * drik.a00
                                         + drij.a10 * drik.a10
                                         + drij.a20 * drik.a20;
                        double cost_ijk = dot_ij_ik / (rij * rik);
                        if (cost_ijk >  1.0) cost_ijk =  1.0;
                        if (cost_ijk < -1.0) cost_ijk = -1.0;

                        // g(cos theta)
                        double u_ang  = d * d + (h + cost_ijk) * (h + cost_ijk);
                        double gik    = gamma * (1.0 + c * c * (1.0 / (d * d) - 1.0 / u_ang));
                        double dg_dcos = gamma * c * c * 2.0 * (h + cost_ijk) / (u_ang * u_ang);

                        // gradients of cos(theta) w.r.t. drij (u) and drik (v)
                        double inv_rij  = 1.0 / rij;
                        double inv_rik  = 1.0 / rik;
                        double inv_rij2 = inv_rij * inv_rij;
                        double inv_rik2 = inv_rik * inv_rik;
                        double inv_rij_rik = inv_rij * inv_rik;

                        Matrix31 term_drij  = drij * (dot_ij_ik * inv_rij2);
                        Matrix31 dcos_ddrij = (drik - term_drij) * inv_rij_rik;

                        Matrix31 term_drik  = drik * (dot_ij_ik * inv_rik2);
                        Matrix31 dcos_ddrik = (drij - term_drik) * inv_rij_rik;

                        // convert to gradients w.r.t. r_i, r_j, r_k (central i)
                        Matrix31 dcos_dri = (dcos_ddrij + dcos_ddrik) * (-1.0);
                        Matrix31 dcos_drj = dcos_ddrij;
                        Matrix31 dcos_drk = dcos_ddrik;

                        Matrix31 dg_dri = dcos_dri * dg_dcos;
                        Matrix31 dg_drj = dcos_drj * dg_dcos;
                        Matrix31 dg_drk = dcos_drk * dg_dcos;

                        // e_ijk = exp(mu * (r_ij - r_ik)), with mu = 2mu from parameter
                        double eijk = std::exp(mu * (rij - rik));

                        // gradients of cutoff and exponential
                        Matrix31 dfc_dri = seik * (-dfc_drik);
                        Matrix31 dfc_drj(0.0, 0.0, 0.0);
                        Matrix31 dfc_drk = seik * dfc_drik;

                        Matrix31 de_dri = (seik - seij) * (mu * eijk);
                        Matrix31 de_drj = seij * (mu * eijk);
                        Matrix31 de_drk = seik * (-mu * eijk);

                        // gradients of X_ijk = fcik * gik * eijk
                        Matrix31 dX_dri = dfc_dri * (gik * eijk)
                                        + dg_dri * (fcik * eijk)
                                        + de_dri * (fcik * gik);
                        Matrix31 dX_drj = dfc_drj * (gik * eijk)
                                        + dg_drj * (fcik * eijk)
                                        + de_drj * (fcik * gik);
                        Matrix31 dX_drk = dfc_drk * (gik * eijk)
                                        + dg_drk * (fcik * eijk)
                                        + de_drk * (fcik * gik);

                        Xij       += fcik * gik * eijk;
                        dXij_dri  = dXij_dri + dX_dri;
                        dXij_drj  = dXij_drj + dX_drj;
                        dXij_drk  = dXij_drk + dX_drk;
                    }
                }

                // ---------- Xji (j-centered) part ----------
                if (rjk > 0.0)
                {
                    // cutoff fc(r_jk)
                    double fcjk, dfc_drjk;
                    if (rjk <= R - D)
                    {
                        fcjk     = 1.0;
                        dfc_drjk = 0.0;
                    }
                    else if (rjk >= R + D)
                    {
                        fcjk     = 0.0;
                        dfc_drjk = 0.0;
                    }
                    else
                    {
                        double xx = 0.5 * M_PI * (rjk - R) / D;
                        fcjk     = 0.5 - 0.5 * std::sin(xx);
                        dfc_drjk = -0.5 * std::cos(xx) * (0.5 * M_PI / D);
                    }

                    if (fcjk > 0.0)
                    {
                        Matrix31 sejk = drjk * (1.0 / rjk);

                        // cos(theta_jik) with central j: angle between r_ji and r_jk
                        double dot_ji_jk = drji.a00 * drjk.a00
                                         + drji.a10 * drjk.a10
                                         + drji.a20 * drjk.a20;
                        double cost_jik = dot_ji_jk / (rji * rjk);
                        if (cost_jik >  1.0) cost_jik =  1.0;
                        if (cost_jik < -1.0) cost_jik = -1.0;

                        double u_ang2   = d * d + (h + cost_jik) * (h + cost_jik);
                        double gjk      = gamma * (1.0 + c * c * (1.0 / (d * d) - 1.0 / u_ang2));
                        double dg2_dcos = gamma * c * c * 2.0 * (h + cost_jik) / (u_ang2 * u_ang2);

                        // gradients of cos(theta) w.r.t. drji (u) and drjk (v)
                        double inv_rji  = 1.0 / rji;
                        double inv_rjk  = 1.0 / rjk;
                        double inv_rji2 = inv_rji * inv_rji;
                        double inv_rjk2 = inv_rjk * inv_rjk;
                        double inv_rji_rjk = inv_rji * inv_rjk;

                        Matrix31 term_drji  = drji * (dot_ji_jk * inv_rji2);
                        Matrix31 dcos_ddrji = (drjk - term_drji) * inv_rji_rjk;

                        Matrix31 term_drjk  = drjk * (dot_ji_jk * inv_rjk2);
                        Matrix31 dcos_ddrjk = (drji - term_drjk) * inv_rji_rjk;

                        // convert to gradients w.r.t. r_j, r_i, r_k (central j)
                        Matrix31 dcos_drj = (dcos_ddrji + dcos_ddrjk) * (-1.0);
                        Matrix31 dcos_dri = dcos_ddrji;
                        Matrix31 dcos_drk = dcos_ddrjk;

                        Matrix31 dg_dri2 = dcos_dri * dg2_dcos;
                        Matrix31 dg_drj2 = dcos_drj * dg2_dcos;
                        Matrix31 dg_drk2 = dcos_drk * dg2_dcos;

                        // e_jik = exp(mu * (r_ji - r_jk)), r_ji = rij
                        double ejik = std::exp(mu * (rji - rjk));

                        // gradients of cutoff and exponential (central j)
                        Matrix31 dfc_dri2(0.0, 0.0, 0.0);
                        Matrix31 dfc_drj2 = sejk * (-dfc_drjk);
                        Matrix31 dfc_drk2 = sejk * dfc_drjk;

                        Matrix31 de_drj2 = (sejk - seji) * (mu * ejik);
                        Matrix31 de_dri2 = seji * (mu * ejik);
                        Matrix31 de_drk2 = sejk * (-mu * ejik);

                        // gradients of X_jik = fcjk * gjk * ejik
                        Matrix31 dX2_dri = dfc_dri2 * (gjk * ejik)
                                         + dg_dri2 * (fcjk * ejik)
                                         + de_dri2 * (fcjk * gjk);
                        Matrix31 dX2_drj = dfc_drj2 * (gjk * ejik)
                                         + dg_drj2 * (fcjk * ejik)
                                         + de_drj2 * (fcjk * gjk);
                        Matrix31 dX2_drk = dfc_drk2 * (gjk * ejik)
                                         + dg_drk2 * (fcjk * ejik)
                                         + de_drk2 * (fcjk * gjk);

                        Xji       += fcjk * gjk * ejik;
                        dXji_dri  = dXji_dri + dX2_dri;
                        dXji_drj  = dXji_drj + dX2_drj;
                        dXji_drk  = dXji_drk + dX2_drk;
                    }
                }

                // store gradients for this k
                dXij_dri_list.push_back(dXij_dri);
                dXij_drj_list.push_back(dXij_drj);
                dXij_drk_list.push_back(dXij_drk);

                dXji_dri_list.push_back(dXji_dri);
                dXji_drj_list.push_back(dXji_drj);
                dXji_drk_list.push_back(dXji_drk);

                k_list.push_back(k);
            } // end k loop

            // bond order terms
            double bij = 1.0 / std::sqrt(1.0 + Xij);
            double bji = 1.0 / std::sqrt(1.0 + Xji);
            double b_sym = 0.5 * (bij + bji);

            // pair contribution to energy
            double u_pair = fcij * (VR - b_sym * VA);

            // split energy equally between i and j
            U_atom[i] += 0.5 * u_pair;
            U_atom[j] += 0.5 * u_pair;

            // --------- forces from bond-order dependence (Xij, Xji) ---------
            // For symmetric ABOP: dU/dXij = fcij * VA * 0.25 / ( (1+Xij)^(3/2) )
            //                     dU/dXji = fcij * VA * 0.25 / ( (1+Xji)^(3/2) )
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

            for (size_t idx = 0; idx < k_list.size(); ++idx)
            {
                int    kk       = k_list[idx];
                Matrix31 dXij_dri = dXij_dri_list[idx];
                Matrix31 dXij_drj = dXij_drj_list[idx];
                Matrix31 dXij_drk = dXij_drk_list[idx];

                Matrix31 dXji_dri = dXji_dri_list[idx];
                Matrix31 dXji_drj = dXji_drj_list[idx];
                Matrix31 dXji_drk = dXji_drk_list[idx];

                // F = - dU/dr = - (dU/dXij * dXij/dr + dU/dXji * dXji/dr)
                Matrix31 Fi_xyz = dXij_dri * (-fk_ij) + dXji_dri * (-fk_ji);
                Matrix31 Fj_xyz = dXij_drj * (-fk_ij) + dXji_drj * (-fk_ji);
                Matrix31 Fk_xyz = dXij_drk * (-fk_ij) + dXji_drk * (-fk_ji);

                data.atoms[i].f = data.atoms[i].f + Fi_xyz;
                data.atoms[j].f = data.atoms[j].f + Fj_xyz;
                data.atoms[kk].f = data.atoms[kk].f + Fk_xyz;
            }

            // --------- pair (radial) forces from fc, VR, VA ---------
            double Fij_scalar = nfcij * (VR - b_sym * VA)
                              + fcij  * (nVR - b_sym * nVA);

            // as in original code: F_i = -dU/dr_i = -Fij * (dr/di) = -Fij*(-seij) = Fij*seij
            Matrix31 fij_xyz = seij * Fij_scalar;

            data.atoms[i].f = data.atoms[i].f + fij_xyz;
            data.atoms[j].f = data.atoms[j].f - fij_xyz;
        }
        
    }

    // total force
    for (i = 0; i < data.n; i++)
    {
        data.F_all = data.F_all + data.atoms[i].f;
    }
}
