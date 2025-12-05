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

void BeW_potential(const parameter1 &pr1, const parameter2 &pr2, Data &data, vector<double> &U_atom)
{

    int i, j, k;
    // int N = data.n;
    // data.U_all = 0.0;
    // double U;
    //  data.K_all = 0.0;
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    //data.f_all = 0.0;

    double D0 = pr2.D0;
    double R = pr2.R;
    double D = pr2.D;
    double mu = pr2.mu; //
    double rf = pr2.rf;
    double bf = pr2.bf;
    double r0 = pr2.r0;
    double beta = pr2.beta;
    double S = pr2.S;
    double gamma = pr2.gamma;
    double c = pr2.c;
    double d = pr2.d;
    double h = pr2.h;

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
            Matrix31 nXij = Matrix31(0.0, 0.0, 0.0);
            Matrix31 Sum_nXij = Matrix31(0.0, 0.0, 0.0);
            Matrix31 nXji = Matrix31(0.0, 0.0, 0.0);
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
            std::vector<Matrix31> nX_list;
            std::vector<Matrix31> seik_list;
            std::vector<Matrix31> sejk_list;
            std::vector<int>      k_list;

            nX_list.clear();
            seik_list.clear();
            sejk_list.clear();
            k_list.clear();

            for (k = 0; k < data.n; k++)
            {
                if (k == i || k == j)
                    continue;
                // E--------------------------------------------------------------------------------------------
                const Atom &ak = data.atoms[k];
                Matrix31 drik = ak.r - ai.r; // distance
                // Matrix31 drjk = ak.r - aj.r;

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
                double rjk2 = (drjk.a00 * drjk.a00 + drjk.a10 * drjk.a10 + drjk.a20 * drjk.a20); // distance

                double rik = sqrt(rik2); // distance
                double rjk = sqrt(rjk2); // distance
                double fcik, fcjk;       // |R-dr| <= D

                double eijk = exp(mu * (rij - rik));
                //double ejik = exp(mu * (rij - rjk));
                // E--------------------------------------------------------------------------------------------
                // F--------------------------------------------------------------------------------------------
                Matrix31 eik = Matrix31(0.0, 1.0, 0.0);
                Matrix31 ejk = Matrix31(0.0, 0.0, 1.0);
                Matrix31 seik = drik * (1.0 / rik);
                Matrix31 sejk = drjk * (1.0 / rjk);

                Matrix31 nfcik, nfcjk;
                // F--------------------------------------------------------------------------------------------
                //------------------------------------------------------------
                //  bij
                //------------------------------------------------------------
                double cost_ijk = (rij2 + rik2 - rjk2) / (2.0 * rij * rik); // theta_ijk angle is i?

                if (cost_ijk > 1.0)
                    cost_ijk = 1.0;
                if (cost_ijk < -1.0)
                    cost_ijk = -1.0;

                double gik = gamma * (1 + c * c * (1.0 / (d * d) - 1 / (d * d + (h + cost_ijk) * (h + cost_ijk))));
                Matrix31 ncost_ijk = (1.0 / (2.0 * rik) - (rik / (2 * rij * rij) + rjk * rjk / (2.0 * rij * rij * rik))) * eij + ((-rij / (2.0 * rik * rik) + 1 / (2.0 * rij) + rjk * rjk / (2.0 * rij * rik * rik))) * eik + ((-rjk / (rij * rik))) * ejk;
                Matrix31 ngik = gamma * c * c * 2 * (h + cost_ijk) * 1.0 / (d * d + (h + cost_ijk) * (h + cost_ijk)) * 1.0 / (d * d + (h + cost_ijk) * (h + cost_ijk)) * ncost_ijk;
                Matrix31 neijk = mu * eijk * (eij - eik);

                if (rik <= R - D) // 1
                {
                    fcik = 1.0;
                    nfcik = Matrix31(0.0, 0.0, 0.0);
                }
                else if (rik >= R + D) // 0
                {
                    fcik = 0.0;
                    nfcik = Matrix31(0.0, 0.0, 0.0);
                }
                else
                {
                    fcik = 0.5 - 0.5 * sin(0.5 * M_PI * (rik - R) / D);
                    nfcik = -1.0 * M_PI * 0.25 * (1.0 / D) * cos(0.5 * M_PI * (rik - R) / D) * eik;
                }
                Xij += fcik * gik * eijk;
                //Sum_Xij = Sum_Xij + Xij;
                nXij = nfcik * (gik * eijk) + fcik * (ngik * eijk + gik * neijk);
                //Sum_nXij = Sum_nXij + nXij;

                nX_list.push_back(nXij);
                seik_list.push_back(seik);
                sejk_list.push_back(sejk);
                k_list.push_back(k);
            }
            double fk =  fcij * VA * (0.5 * 1.0 / (sqrt(1 + Xij) * (1 + Xij)));
            for(int k1 = 0; k1 < nX_list.size(); k1++)
            {
                int k = k_list[k1];
                Matrix31 nXij = nX_list[k1];
                Matrix31 seik = seik_list[k1];
                Matrix31 sejk = sejk_list[k1];

                Matrix31 fk_ijk =  fk * nXij;
                //f. ijk to xyz transform
                //-gradU. -->  -(du/dik * dik/dk) - (du/djk * djk/dk) 
                Matrix31 fk_xyz = - fk_ijk.a10 * seik - fk_ijk.a20 * sejk;//fk. ijk to xyz transform
                //-gradU. -->  -(du/dij * dij/di) - (du/dik * dik/di)
                Matrix31 fi_xyz = - fk_ijk.a00 * (-1.0 * seij) - fk_ijk.a10 *(-1.0 * seik);//atom_k to atom_i force
                //-gradU. -->  -(du/dij * dij/dj) - (du/djk * djk/dj)
                Matrix31 fj_xyz = - fk_ijk.a00 * seij - fk_ijk.a20 * (-1.0 * sejk);//atom_k to atom_j force  && d/di = -d/dj
                data.atoms[k].f = data.atoms[k].f + 0.5*fk_xyz;
                data.atoms[i].f = data.atoms[i].f + 0.5*fi_xyz;
                data.atoms[j].f = data.atoms[j].f + 0.5*fj_xyz;

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
