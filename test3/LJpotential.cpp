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
void LJ_potential(Data &data, const parameter1 &p1, vector<double> &U_atom,
                  Matrix31 &F)
{
    int i, j;
    int N = data.n;
    data.U_all = 0.0;
    data.K_all = 0.0;
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    double epsilon = p1.epsilon;
    double sigma = p1.sigma;
    double m = p1.m;


    U_atom.assign(N, 0.0);
    data.U_all = 0.0;
    F = Matrix31(0.0, 0.0, 0.0);

    for (i = 0; i < data.n; i++)
    {
        const Atom &ai = data.atoms[i];
        data.atoms[i].dv = Matrix31(0.0, 0.0, 0.0);

        for (j = i + 1; j < data.n; j++)
        {
            const Atom &aj = data.atoms[j];
            
            Matrix31 dr = aj.r - ai.r;
            // double dx = aj.r.a00 - ai.r.a00;// after translate to matrix
            // double dy = aj.r.a10 - ai.r.a10;
            // double dz = aj.r.a20 - ai.r.a20;

            double r = sqrt(dr.a00 * dr.a00 + dr.a10 * dr.a10 + dr.a20 * dr.a20);

            // fix \\

            // double r2 = r * r;
            // double r6 = r2 * r2 * r2;
            // double r12 = r6 * r6;

            // double uij = 4.0 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));//pow(a,b) = a^b
            // data.U_all += uij;

            U_atom[i] += 0.5 * uij;
            U_atom[j] += 0.5 * uij;
            //dU/dr
            F = 24 * epsilon *(2* pow(sigma / r, 12) - pow(sigma / r, 6)) * (dr * (1/(r*r)));\

            data.F_all = F + data.F_all;
            //F = coeff * dr + F;

            data.atoms[i].dv = F * (1/m) + data.atoms[i].dv;
            data.atoms[j].dv = -1 * F * (1/m) + data.atoms[j].dv;
        };
    };
}