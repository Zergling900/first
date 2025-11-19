#include <iostream>
#include <cmath>
#include <vector>

#include "3.h"
#include "void.h"

using namespace std;

// ----------------------------
// time
// ----------------------------
void evolution(const parameter1 &pr1, Data &Data0,vector<double> &U_atom)
{
    int N = Data0.n;
    double dt = pr1.dt;
    double hdt = dt / 2.0;
    double epsilon = pr1.epsilon;
    double kb = pr1.kb;
    double m1 = 1.0 / pr1.m;
    double T = pr1.T;
    double sigma = pr1.sigma;
    double steps = pr1.steps;

    //vector<double> U_atom(N);

    Matrix31 r0 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 p0 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 f0 = Matrix31(0.0, 0.0, 0.0);
    // Matrix31 r1 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 p1 = Matrix31(0.0, 0.0, 0.0);
    // Matrix31 f1 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 r2 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 p2 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 f2 = Matrix31(0.0, 0.0, 0.0);

    for (int i = 0; i < N; i++)
    {
        // 0 is (t)
        // 1 is (t+hdt)
        // 2 is (t+dt)
        r0 = Data0.atoms[i].r;
        p0 = Data0.atoms[i].p;
        f0 = Data0.atoms[i].f;

        p1 = p0 + hdt * f0;
        r2 = r0 + dt * p1 * m1;

        if (r2.a00 > Data0.Box.a00)
            r2.a00 = r2.a00 - Data0.Box.a00;
        if (r2.a00 < 0.0)
            r2.a00 = r2.a00 + Data0.Box.a00;
        if (r2.a10 > Data0.Box.a11)
            r2.a10 = r2.a10 - Data0.Box.a11;
        if (r2.a10 < 0.0)
            r2.a10 = r2.a10 + Data0.Box.a11;
        if (r2.a20 > Data0.Box.a22)
            r2.a20 = r2.a20 - Data0.Box.a22;
        if (r2.a20 < 0.0)
            r2.a20 = r2.a20 + Data0.Box.a22;

        Data0.atoms[i].r = r2;
        Data0.atoms[i].p = p1;
    }

    LJ_potential(Data0, pr1, U_atom);

    for (int i = 0; i < N; i++)
    {
        p1 = Data0.atoms[i].p;
        f2 = Data0.atoms[i].f;
        p2 = p1 + dt * f2 * 0.5;
        Data0.atoms[i].p = p2;
    }
    energy(Data0, pr1, U_atom);
}