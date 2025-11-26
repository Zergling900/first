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
    //int N = data.n;
    //data.U_all = 0.0;
    //data.K_all = 0.0;
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    data.f_all = 0.0;
    double epsilon = pr1.epsilon * pr1.kb;
    double sigma = pr1.sigma;
    double sigma2 = sigma * sigma;

    U_atom.assign(data.n, 0.0);//initialize

    Matrix31 F = Matrix31(0.0, 0.0, 0.0);

    //PBC,but this just for positive crystal (test)
    // double Lx = sqrt(data.Box.a00 * data.Box.a00 + data.Box.a01 * data.Box.a01 + data.Box.a02 * data.Box.a02);
    // double Ly = sqrt(data.Box.a10 * data.Box.a10 + data.Box.a11 * data.Box.a11 + data.Box.a12 * data.Box.a12);
    // double Lz = sqrt(data.Box.a20 * data.Box.a20 + data.Box.a21 * data.Box.a21 + data.Box.a22 * data.Box.a22);
    // double Lxh = Lx / 2.0;
    // double Lyh = Ly / 2.0;
    // double Lzh = Lz / 2.0;
     double Lx = data.Box.a00;
     double Ly = data.Box.a11;
     double Lz = data.Box.a22;
     double Lxh = data.Box.a00 / 2.0;
     double Lyh = data.Box.a11 / 2.0;
     double Lzh = data.Box.a22 / 2.0;

    double rc2 = 2.5 * 2.5 * sigma2;//rc^2 , p30,  2.5~3.5 * sigma

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

            //if < or > L/2 , translate to -L/2
            //------------------------------------------------------------------ 

            if (dr.a00 >  Lxh) dr.a00 -= Lx;
            if (dr.a00 < -Lxh) dr.a00 += Lx;

            if (dr.a10 >  Lyh) dr.a10 -= Ly;
            if (dr.a10 < -Lyh) dr.a10 += Ly;

            if (dr.a20 >  Lzh) dr.a20 -= Lz;
            if (dr.a20 < -Lzh) dr.a20 += Lz;
            //------------------------------------------------------------------

            // this r is |dr|!
            //double r = sqrt(dr.a00 * dr.a00 + dr.a10 * dr.a10 + dr.a20 * dr.a20);
            double r2 = dr.a00 * dr.a00 + dr.a10 * dr.a10 + dr.a20 * dr.a20;
            if (r2 == 0.0)
            {
                continue;
                //r2 = 1e-14;
            }

            if (r2 > rc2) continue;// if r > rc, skip

            double sr2 = sigma2 / r2;// sigma^2/r^2
            double sr6 = sr2 * sr2 * sr2;
            double sr12 = sr6 * sr6;

            double uij = 4 * epsilon * (sr12 - sr6);//epsilon is epsilon/kb
            //data.U_all += uij;//in energy

            U_atom[i] += 0.5 * uij;
            U_atom[j] += 0.5 * uij;

            // dU/dr,this r is |dr|,so need /r^2
            // F is 3*1 vector
            F = (-48.0) * epsilon * (sr12 - 0.5 * sr6) / (r2) * dr;

            //data.F_all = F + data.F_all;//mistake
            //if slow, can delete this part,or change to f^2

            data.atoms[i].f = data.atoms[i].f + F;
            data.atoms[j].f = data.atoms[j].f - F;
        };
    };

    for (i = 0; i < data.n; i++)
    {
        data.F_all = data.F_all + data.atoms[i].f;
    };
}