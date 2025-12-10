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
            // PBC*******************************************

            double rij2 = (drij.a00 * drij.a00 + drij.a10 * drij.a10 + drij.a20 * drij.a20); // distance
            double rij = sqrt(rij2);                                                         // distance
            double bsR = beta * sqrt(2 * S);
            double bsA = beta * sqrt(2 / S);
            
            double VR = (D0 / (S - 1.0)) * exp(-bsR * (rij - r0));
            double VA = (S * D0 / (S - 1.0)) * exp(-bsA * (rij - r0));
            
            // double fc1 = 1.0;                                          // r <= R-D
            double fcij; // |R-dr| <= D
            // double fc3 = 0.0;                                          // r >= R+D
            double Xij = 0.0, Xji = 0.0;
            //double Sum_Xij = 0.0;
            // E--------------------------------------------------------------------------------------------
            // F--------------------------------------------------------------------------------------------
            Matrix31 seij = drij * (1.0 / rij);
            //Matrix31 eij = Matrix31(1.0, 0.0, 0.0);
            double nfcij;
            double nVR = D0 / (S - 1.0) * (-bsR) * exp(-bsR * (rij - r0));
            double nVA = S * D0 / (S - 1.0) * (-bsA) * exp(-bsA * (rij - r0));
            // F--------------------------------------------------------------------------------------------

            //fc_ij------------------------------------------------
            if (rij <= (R - D)) // 1
            {
                fcij = 1.0;
                nfcij = 0.0;
            }
            else if (rij >= (R + D)) // 0
            {
                fcij = 0.0;
                nfcij = 0.0;
            }
            else // 1/2 -1/2sin(pi/2(r-R)/D)
            {
                fcij = 0.5 - 0.5 * sin(0.5 * M_PI * (rij - R) / D);
                nfcij = -1.0 * M_PI * 0.25 * (1.0 / D) * cos(0.5 * M_PI * (rij - R) / D);
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
                else if (drik.a00 >= Lxh)
                    drik.a00 -= Lx;
                if (drik.a10 < -Lyh)
                    drik.a10 += Ly;
                else if (drik.a10 >= Lyh)
                    drik.a10 -= Ly;
                if (drik.a20 < -Lzh)
                    drik.a20 += Lz;
                else if (drik.a20 >= Lzh)
                    drik.a20 -= Lz;
                // PBC*******************************************

                //Matrix31 drjk = drik - drij;

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

                //dcos = drik/(rij*rik) - (drij \cdot drik)/(rij^3*rik) * drij

                Matrix31 term_drij = drij * (dot_ij_ik * inv_rij2);
                Matrix31 dcos_ddrij = (drik - term_drij) * inv_rij_rik;

                Matrix31 term_drik = drik * (dot_ij_ik * inv_rik2);
                Matrix31 dcos_ddrik = (drij - term_drik) * inv_rij_rik;

                // convert to gradients w.r.t. r_i, r_j, r_k
                //dcos = dcos/drij *drij/dri + dcos/drik *drik/dri
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
            //!ZBL?
            Uij = Uij + fcij * (VR - (bij)*VA);
            //nUij = nUij + nfcij * (VR - (bij)*VA) + fcij * (nVR - (nbij)*VA - (bij)*nVA);
            double Fij = nfcij * ( VR - bij * VA) + fcij * (nVR - bij * nVA);
            //gradU. -->  -(du/dij * dij/di) 
            Matrix31 fij_xyz = - Fij * (-1.0 * seij);//atom_k to atom_i force
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
        U_atom[i] += 0.5 * Uij;
        //data.atoms[i].f = nUij;
        // data.atoms[i].f = -1.0 * nUij; // F = -grad U
    } // i

    for (i = 0; i < data.n; i++)
    {
        data.F_all = data.F_all + data.atoms[i].f;
    };
}

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
            // double rf  = p->rf; // ZBL-related, not used yet
            // double bf  = p->bf;

            // Minimum-image displacement r_ij = r_j - r_i
            Matrix31 drij = aj.r - ai.r;

            if (drij.a00 < -Lxh)      drij.a00 += Lx;
            else if (drij.a00 >= Lxh) drij.a00 -= Lx;

            if (drij.a10 < -Lyh)      drij.a10 += Ly;
            else if (drij.a10 >= Lyh) drij.a10 -= Ly;

            if (drij.a20 < -Lzh)      drij.a20 += Lz;
            else if (drij.a20 >= Lzh) drij.a20 -= Lz;

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

            // Radial derivatives of VR, VA
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
                if (drik.a00 < -Lxh)      drik.a00 += Lx;
                else if (drik.a00 >= Lxh) drik.a00 -= Lx;

                if (drik.a10 < -Lyh)      drik.a10 += Ly;
                else if (drik.a10 >= Lyh) drik.a10 -= Ly;

                if (drik.a20 < -Lzh)      drik.a20 += Lz;
                else if (drik.a20 >= Lzh) drik.a20 -= Lz;

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

                // ---------- Xij (i-centered) part ----------
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
