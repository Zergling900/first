#include <iostream>
#include <cmath>
#include <vector>

#include "3.h"
#include "void.h"

using namespace std;

// ----------------------------
// time
// ----------------------------
void BeW_evolution1(const parameter1 &pr1, const parameter2 &pr2_WW,const parameter2 &pr2_BB,const parameter2 &pr2_WB, Data &data,vector<double> &U_atom)
{
    int N = data.n;
    double dt = pr1.dt;
    double hdt = dt / 2.0;

    double kb = pr1.kb;
    double Q = pr1.Q;
    int g = pr1.g;
    double s0 = data.s0;
    double ps0 = data.ps0;
    double T = pr1.T;
    double H0 = pr1.H0;
    
    double mw = pr1.mw;
    double mb = pr1.mb;
    //vector<double> U_atom(N);

    Matrix31 r0 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 p0 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 f0 = Matrix31(0.0, 0.0, 0.0);
    double u0 = 0.0, u1 = 0.0, u2 = 0.0;

    Matrix31 r1 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 p1 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 f1 = Matrix31(0.0, 0.0, 0.0);

    Matrix31 r2 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 p2 = Matrix31(0.0, 0.0, 0.0);
    Matrix31 f2 = Matrix31(0.0, 0.0, 0.0);

    double s1 = 0.0,s2 = 0.0, ps1 = 0.0, ps2 = 0.0, ps3 = 0.0, ps4 = 0.0, ps5 = 0.0;
    

    //***a***
    s1 = s0 * (1.0 + hdt * ps0 / (2.0 * Q))*(1.0 + hdt * ps0 / (2.0 * Q));
    ps1 = ps0 / (1 +  hdt * ps0 / (2.0 * Q));

    //***b***
    u0 = data.U_all;
    ps2 = ps1 - hdt * u0;
    for (int i = 0; i < N; i++)
    {
        // 0 is (t)
        // 1 is (t+hdt)
        // 2 is (t+dt)
        r0 = data.atoms[i].r;
        p0 = data.atoms[i].p;
        f0 = data.atoms[i].f;

        //!p1 = p0 + hdt * s0 *(f0);
        p1 = p0 + hdt * s1 *(f0);

    //***c***
        if (data.atoms[i].name == "W")
        {
            r2 = r0 + (dt * p1 * (1.0/ (mw * s1)));
        }
        else if(data.atoms[i].name == "Be")
        {
            r2 = r0 + (dt * p1 * (1.0/ (mb * s1)));
        }

        if (r2.a00 >= data.Box.a00)
            r2.a00 = r2.a00 - data.Box.a00;
        if (r2.a00 < 0.0)
            r2.a00 = r2.a00 + data.Box.a00;
        if (r2.a10 >= data.Box.a11)
            r2.a10 = r2.a10 - data.Box.a11;
        if (r2.a10 < 0.0)
            r2.a10 = r2.a10 + data.Box.a11;
        if (r2.a20 >= data.Box.a22)
            r2.a20 = r2.a20 - data.Box.a22;
        if (r2.a20 < 0.0)
            r2.a20 = r2.a20 + data.Box.a22;

        data.atoms[i].r = r2;
        data.atoms[i].p = p1;
    }

    double K = 0.0;
    for (int i = 0; i < data.n; i++)
    {
        if(data.atoms[i].name == "W")
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 /(pr1.mw * s1 * s1));
        else if(data.atoms[i].name == "Be")
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 / (pr1.mb * s1 * s1));
    }
    ps3 = ps2 - dt *(-K + g *kb * T * log(s1) - H0 + g *kb * T);
    //***d***
    //---
    BeW_potential2(pr1, pr2_WW,pr2_BB,pr2_WB, data, U_atom);
    
    for (int i = 0; i < N; i++)
    {
        p1 = data.atoms[i].p;
        f2 = data.atoms[i].f;
        p2 = p1 + hdt * s1 *(f2);
        data.atoms[i].p = p2;
    }
    energy(data, pr1, U_atom);
    u2 = data.U_all;
    ps4 = ps3 - hdt * u2;
    //***e***
    s2 = s1 * (1.0 + hdt * ps4 / (2.0 * Q))*(1.0 + hdt * ps4 / (2.0 * Q));
    ps5 = ps4 /(1.0+ hdt * ps4 / (2.0 * Q));


    //***new***
    data.s0 = s2;
    data.ps0 = ps5;
    //***new***

}
