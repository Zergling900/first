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
void LJ_potential(Data &data, const parameter1 &pr1, vector<double> &U_atom)
{
    int i, j;
    // int N = data.n;
    // data.U_all = 0.0;
    // data.K_all = 0.0;
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    //data.f_all = 0.0;
    double epsilon = pr1.epsilon * pr1.kb;
    double sigma = pr1.sigma;
    double sigma2 = sigma * sigma;

    U_atom.assign(data.n, 0.0); // initialize

    Matrix31 F = Matrix31(0.0, 0.0, 0.0);

    // PBC,but this just for positive crystal (test)
    //  double Lx = sqrt(data.Box.a00 * data.Box.a00 + data.Box.a01 * data.Box.a01 + data.Box.a02 * data.Box.a02);
    //  double Ly = sqrt(data.Box.a10 * data.Box.a10 + data.Box.a11 * data.Box.a11 + data.Box.a12 * data.Box.a12);
    //  double Lz = sqrt(data.Box.a20 * data.Box.a20 + data.Box.a21 * data.Box.a21 + data.Box.a22 * data.Box.a22);
    //  double Lxh = Lx / 2.0;
    //  double Lyh = Ly / 2.0;
    //  double Lzh = Lz / 2.0;
    double Lx = data.Box.a00;
    double Ly = data.Box.a11;
    double Lz = data.Box.a22;
    double Lxh = data.Box.a00 / 2.0;
    double Lyh = data.Box.a11 / 2.0;
    double Lzh = data.Box.a22 / 2.0;

    double rc2 = 2.5 * 2.5 * sigma2; // rc^2 , p30,  2.5~3.5 * sigma

    for (i = 0; i < data.n; i++)
    {
        data.atoms[i].f = Matrix31(0.0, 0.0, 0.0);
    };

    for (i = 0; i < data.n; i++)
    {
        const Atom &ai = data.atoms[i];

        for (j = i + 1; j < data.n; j++)
        {
            const Atom &aj = data.atoms[j];

            Matrix31 dr = aj.r - ai.r;
            // double dx = aj.r.a00 - ai.r.a00;// after translate to matrix
            // double dy = aj.r.a10 - ai.r.a10;
            // double dz = aj.r.a20 - ai.r.a20;

            // if < or > L/2 , translate to -L/2
            //------------------------------------------------------------------

            if (dr.a00 > Lxh)
                dr.a00 -= Lx;
            if (dr.a00 < -Lxh)
                dr.a00 += Lx;

            if (dr.a10 > Lyh)
                dr.a10 -= Ly;
            if (dr.a10 < -Lyh)
                dr.a10 += Ly;

            if (dr.a20 > Lzh)
                dr.a20 -= Lz;
            if (dr.a20 < -Lzh)
                dr.a20 += Lz;
            //------------------------------------------------------------------

            // this r is |dr|!
            // double r = sqrt(dr.a00 * dr.a00 + dr.a10 * dr.a10 + dr.a20 * dr.a20);
            double r2 = dr.a00 * dr.a00 + dr.a10 * dr.a10 + dr.a20 * dr.a20;
            if (r2 == 0.0)
            {
                continue;
                // r2 = 1e-14;
            }

            if (r2 > rc2)
                continue; // if r > rc, skip

            double sr2 = sigma2 / r2; // sigma^2/r^2
            double sr6 = sr2 * sr2 * sr2;
            double sr12 = sr6 * sr6;

            double uij = 4 * epsilon * (sr12 - sr6); // epsilon is epsilon/kb
            // data.U_all += uij;//in energy

            U_atom[i] += 0.5 * uij;
            U_atom[j] += 0.5 * uij;

            // dU/dr,this r is |dr|,so need /r^2
            // F is 3*1 vector
            F = (-48.0) * epsilon * (sr12 - 0.5 * sr6) / (r2)*dr;

            // data.F_all = F + data.F_all;//mistake
            // if slow, can delete this part,or change to f^2

            data.atoms[i].f = data.atoms[i].f + F;
            data.atoms[j].f = data.atoms[j].f - F;
        };
    };

    for (i = 0; i < data.n; i++)
    {
        data.F_all = data.F_all + data.atoms[i].f;
    };
}

void BeW_potential(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,
                   vector<double> &U_atom)
{

    int i, j, k;
    // int N = data.n;
    // data.U_all = 0.0;
    // double U;
    //  data.K_all = 0.0;
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    //data.f_all = 0.0;

    // double m = pr1.m;

    U_atom.assign(data.n, 0.0); // initialize

    // Matrix31 F = Matrix31(0.0, 0.0, 0.0);
    //--------------------------------------------------------
    double Lx = data.Box.a00;
    double Ly = data.Box.a11;
    double Lz = data.Box.a22;
    double Lxh = data.Box.a00 / 2.0;
    double Lyh = data.Box.a11 / 2.0;
    double Lzh = data.Box.a22 / 2.0;

    for (i = 0; i < data.n; i++)
    {
        data.atoms[i].f = Matrix31(0.0, 0.0, 0.0);
    };

    for (i = 0; i < data.n; i++)
    {
        const Atom &ai = data.atoms[i];

        double Uij = 0.0;
        Matrix31 nUij = Matrix31(0.0, 0.0, 0.0);

        for (j = 0; j < data.n; j++)
        {
            if (i == j)
                continue;
            // E--------------------------------------------------------------------------------------------
            const Atom &aj = data.atoms[j];
            Matrix31 drij = aj.r - ai.r; // distance
            // Choose ABOP parameter set based on species of i and j
            const parameter2 *p = nullptr;
            const std::string &ni = ai.name;
            const std::string &nj = aj.name;

            // Here we assume:
            //  - pr2_WW: parameters for W-W pairs
            //  - pr2_BB: parameters for Be-Be pairs
            //  - pr2_WB: parameters for W-Be (or Be-W) pairs
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
            double D0   = p->D0;
            double R    = p->R;
            double D    = p->D;
            double mu   = p->mu;
            double rf   = p->rf;
            double bf   = p->bf;
            double r0   = p->r0;
            double beta = p->beta;
            double S    = p->S;
            double gamma = p->gamma;
            double c     = p->c;
            double d     = p->d;
            double h     = p->h;
            // PBC*******************************************
            if (drij.a00 < -Lxh)
                drij.a00 += Lx;
            else if (drij.a00 > Lxh)
                drij.a00 -= Lx;
            if (drij.a10 < -Lyh)
                drij.a10 += Ly;
            else if (drij.a10 > Lyh)
                drij.a10 -= Ly;
            if (drij.a20 < -Lzh)
                drij.a20 += Lz;
            else if (drij.a20 > Lzh)
                drij.a20 -= Lz;
            // PBC*******************************************

            double rij2 = (drij.a00 * drij.a00 + drij.a10 * drij.a10 + drij.a20 * drij.a20); // distance
            double rij = sqrt(rij2);                                                         // distance
            double bsR = beta * sqrt(2 * S);
            double bsA = beta * sqrt(2 / S);
            double VR = D0 / (S - 1.0) * exp(-bsR * (rij - r0));
            double VA = S * D0 / (S - 1.0) * exp(-bsA * (rij - r0));

            // double fc1 = 1.0;                                          // r <= R-D
            double fcij; // |R-dr| <= D
            // double fc3 = 0.0;                                          // r >= R+D
            double Xij = 0.0, Xji = 0.0;
            //double Sum_Xij = 0.0;
            // E--------------------------------------------------------------------------------------------
            // F--------------------------------------------------------------------------------------------
            Matrix31 seij = drij * (1.0 / rij);
            Matrix31 eij = Matrix31(1.0, 0.0, 0.0);
            Matrix31 nfcij;
            Matrix31 nVR = D0 / (S - 1.0) * (-bsR) * exp(-bsR * (rij - r0)) * eij;
            Matrix31 nVA = S * D0 / (S - 1.0) * (-bsA) * exp(-bsA * (rij - r0)) * eij;
            // F--------------------------------------------------------------------------------------------

            //fc_ij------------------------------------------------
            if (rij <= R - D) // 1
            {
                fcij = 1.0;
                nfcij = Matrix31(0.0, 0.0, 0.0);
            }
            else if (rij >= R + D) // 0
            {
                fcij = 0.0;
                nfcij = Matrix31(0.0, 0.0, 0.0);
            }
            else // 1/2 -1/2sin(pi/2(r-R)/D)
            {
                fcij = 0.5 - 0.5 * sin(0.5 * M_PI * (rij - R) / D);
                nfcij = -1.0 * M_PI * 0.25 * (1.0 / D) * cos(0.5 * M_PI * (rij - R) / D) * eij;
            }
            //fc_ij------------------------------------------------

            // k
            std::vector<Matrix31> dX_dri_list;
            std::vector<Matrix31> dX_drj_list;
            std::vector<Matrix31> dX_drk_list;
            std::vector<int>      k_list;

            dX_dri_list.clear();
            dX_drj_list.clear();
            dX_drk_list.clear();
            k_list.clear();

            for (k = 0; k < data.n; k++)
            {
                if (k == i || k == j)
                    continue;
                // E--------------------------------------------------------------------------------------------
                const Atom &ak = data.atoms[k];
                Matrix31 drik = ak.r - ai.r; // distance

                // PBC*******************************************
                if (drik.a00 < -Lxh)
                    drik.a00 += Lx;
                else if (drik.a00 > Lxh)
                    drik.a00 -= Lx;
                if (drik.a10 < -Lyh)
                    drik.a10 += Ly;
                else if (drik.a10 > Lyh)
                    drik.a10 -= Ly;
                if (drik.a20 < -Lzh)
                    drik.a20 += Lz;
                else if (drik.a20 > Lzh)
                    drik.a20 -= Lz;
                // PBC*******************************************

                Matrix31 drjk = drik - drij;

                double rik2 = (drik.a00 * drik.a00 + drik.a10 * drik.a10 + drik.a20 * drik.a20); // distance
                double rik = sqrt(rik2); // distance
                if (rik == 0.0)
                {
                    continue;
                }

                // cutoff fc(r_ik)
                double fcik;
                double dfc_drik; // d fcik / d r_ik

                if (rik <= R - D) // 1
                {
                    fcik = 1.0;
                    dfc_drik = 0.0;
                }
                else if (rik >= R + D) // 0
                {
                    fcik = 0.0;
                    dfc_drik = 0.0;
                }
                else
                {
                    fcik = 0.5 - 0.5 * sin(0.5 * M_PI * (rik - R) / D);
                    dfc_drik = -1.0 * M_PI * 0.25 * (1.0 / D) * cos(0.5 * M_PI * (rik - R) / D);
                }

                // unit vectors
                Matrix31 seik = drik * (1.0 / rik);

                // cos(theta_ijk) using vector formula
                double dot_ij_ik = drij.a00 * drik.a00 + drij.a10 * drik.a10 + drij.a20 * drik.a20;
                double cost_ijk = dot_ij_ik / (rij * rik);

                if (cost_ijk > 1.0)
                    cost_ijk = 1.0;
                if (cost_ijk < -1.0)
                    cost_ijk = -1.0;

                // g_ik(costheta)
                double gik = gamma * (1 + c * c * (1.0 / (d * d) - 1.0 / (d * d + (h + cost_ijk) * (h + cost_ijk))));

                // dg/dcos
                double u = d * d + (h + cost_ijk) * (h + cost_ijk);
                double dg_dcos = gamma * c * c * 2.0 * (h + cost_ijk) / (u * u);

                // gradients of cos(theta_ijk) w.r.t. u = drij and v = drik
                double inv_rij = 1.0 / rij;
                double inv_rik = 1.0 / rik;
                double inv_rij2 = inv_rij * inv_rij;
                double inv_rik2 = inv_rik * inv_rik;
                double inv_rij_rik = inv_rij * inv_rik;

                //Matrix31 u = drij;
                //Matrix31 v = drik;

                //double uv_dot = dot_ij_ik;

                Matrix31 term_drij = drij * (dot_ij_ik * inv_rij2);
                Matrix31 dcos_ddrij = (drik - term_drij) * inv_rij_rik;

                Matrix31 term_drik = drik * (dot_ij_ik * inv_rik2);
                Matrix31 dcos_ddrik = (drij - term_drik) * inv_rij_rik;

                // convert to gradients w.r.t. r_i, r_j, r_k
                Matrix31 dcos_dri = (dcos_ddrij + dcos_ddrik) * (-1.0);
                Matrix31 dcos_drj = dcos_ddrij;
                Matrix31 dcos_drk = dcos_ddrik;

                Matrix31 dg_dri = dcos_dri * dg_dcos;
                Matrix31 dg_drj = dcos_drj * dg_dcos;
                Matrix31 dg_drk = dcos_drk * dg_dcos;

                // e_ijk = exp(mu * (r_ij - r_ik))
                double eijk = exp(mu * (rij - rik));

                // gradients of cutoff and exponential
                Matrix31 dfc_dri = seik * (-dfc_drik);
                Matrix31 dfc_drj = Matrix31(0.0, 0.0, 0.0);
                Matrix31 dfc_drk = seik * dfc_drik;

                Matrix31 de_dri = (seik - seij) * (mu * eijk);
                Matrix31 de_drj = seij * (mu * eijk);
                Matrix31 de_drk = seik * (-mu * eijk);

                // gradients of X_ijk = fcik * gik * eijk
                Matrix31 dX_dri = dfc_dri * (gik * eijk) + dg_dri * (fcik * eijk) + de_dri * (fcik * gik);
                Matrix31 dX_drj = dfc_drj * (gik * eijk) + dg_drj * (fcik * eijk) + de_drj * (fcik * gik);
                Matrix31 dX_drk = dfc_drk * (gik * eijk) + dg_drk * (fcik * eijk) + de_drk * (fcik * gik);

                Xij += fcik * gik * eijk;

                dX_dri_list.push_back(dX_dri);
                dX_drj_list.push_back(dX_drj);
                dX_drk_list.push_back(dX_drk);
                k_list.push_back(k);
            }
            double fk = fcij * VA * (0.5 * 1.0 / (sqrt(1 + Xij) * (1 + Xij)));
            for (int k1 = 0; k1 < static_cast<int>(k_list.size()); k1++)
            {
                int k = k_list[k1];

                // forces from X_ij through X_ijk contributions
                Matrix31 dX_dri = dX_dri_list[k1];
                Matrix31 dX_drj = dX_drj_list[k1];
                Matrix31 dX_drk = dX_drk_list[k1];

                Matrix31 Fi_xyz = dX_dri * (-fk); // F = -dU/dr
                Matrix31 Fj_xyz = dX_drj * (-fk);
                Matrix31 Fk_xyz = dX_drk * (-fk);

                // 0.5 factor to avoid double counting i-j pairs
                data.atoms[i].f = data.atoms[i].f + 0.5 * Fi_xyz;
                data.atoms[j].f = data.atoms[j].f + 0.5 * Fj_xyz;
                data.atoms[k].f = data.atoms[k].f + 0.5 * Fk_xyz;
            }
            // k

            double bij = 1.0 / sqrt(1 + Xij);
            /*
            double bji = 1.0 / sqrt(1 + Xji);
            */

            // double bij = 1.0 ;
            // double bji = 1.0 ;

            //Matrix31 nbij = -0.5 * 1.0 / (sqrt(1 + Xij) * (1 + Xij)) * Sum_nXij;

            Uij = Uij + fcij * (VR - (bij)*VA);
            //nUij = nUij + nfcij * (VR - (bij)*VA) + fcij * (nVR - (nbij)*VA - (bij)*nVA);
            Matrix31 Fij = nfcij * ( VR - bij * VA) + fcij * (nVR - bij * nVA);
            //gradU. -->  -(du/dij * dij/di) 
            Matrix31 fij_xyz = - Fij.a00 * (-1.0 * seij);//atom_k to atom_i force
            //-gradU. -->  (du/dij * dij/di) 
            //Matrix31 fj_xyz = - Fij.a00 * seij;//atom_k to atom_j force  && d/di = -d/dj
            data.atoms[i].f = data.atoms[i].f + 0.5*fij_xyz;
            data.atoms[j].f = data.atoms[j].f - 0.5*fij_xyz;
            //* for i != j so * 0.5 

        } // j-------------------------------------------------------------

        // U_atom[i] += 0.5 * Uij;
        // U_atom[j] += 0.5 * Uij;

        // data.atoms[i].f = data.atoms[i].f + nUij;
        // data.atoms[j].f = data.atoms[j].f - 1.0 * nUij;
        U_atom[i] = 0.5 * Uij;
        //data.atoms[i].f = nUij;
        // data.atoms[i].f = -1.0 * nUij; // F = -grad U
    } // i

    for (i = 0; i < data.n; i++)
    {
        data.F_all = data.F_all + data.atoms[i].f;
    };
}
