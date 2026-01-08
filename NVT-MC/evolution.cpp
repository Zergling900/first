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
    //double epsilon = pr1.epsilon;
    //double kb = pr1.kb;
    double mw1 = 1.0 / pr1.mw;
    double mb1 = 1.0 / pr1.mb;
    //double T = pr1.T;
    //double sigma = pr1.sigma;
    //double steps = pr1.steps;
    //!*****new******
    double kb = pr1.kb;
    double Q = pr1.Q/(pr1.E0 * pr1.tau * pr1.tau);
    int g = pr1.g;
    double s0 = data.s0, ps0 = data.ps0;
    // if (s0 == pr1.s0 && ps0 == pr1.ps0)
    //     s0 = pr1.s0;
    //     ps0 = pr1.ps0;
    double T = pr1.T;
    //double dxi = 2.0 / Q *(K - 3.0/2.0 * N * kb * T);
    
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
    //***new***
    double s1, s2, ps2;

    //***new***
    s1 = s0 + hdt * ps0 / Q;
    
    for (int i = 0; i < N; i++)
    {
        // 0 is (t)
        // 1 is (t+hdt)
        // 2 is (t+dt)
        r0 = data.atoms[i].r;
        p0 = data.atoms[i].p;
        f0 = data.atoms[i].f;
        
        //***new***
        //dp/dt = f - (ps/Q)*p
        if (data.atoms[i].name == "W")
        {
            p1 = p0 + hdt * (f0 - (ps0/Q)*p0);
            //s1 = s0 + hdt * ps0 / Q;
            r2 = r0 + (dt * p1 * mw1 ) * (1.0/ (s1 * s1));
            //ps2 = ps0 + dt *(2.0*K/(s1*s1) - g * kb * T / s1);
        }
        else if(data.atoms[i].name == "Be")
        {
            p1 = p0 + hdt * (f0 - (ps0/Q)*p0);
            //s1 = s0 + hdt * ps0 / Q;
            r2 = r0 + (dt * p1 * mb1 ) * (1.0/ (s1 * s1));
            //ps2 = ps0 + dt *(2.0*K/(s1*s1) - g * kb * T / s1);
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
    //***new***
    double K = 0.0;
    for (int i = 0; i < data.n; i++)
    {
        if(data.atoms[i].name == "W")
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 / pr1.mw);
        else if(data.atoms[i].name == "Be")
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                        + data.atoms[i].p.a10 * data.atoms[i].p.a10
                        + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                        * (0.5 / pr1.mb);
    }
    ps2 = ps0 + dt *(2.0*K/(s1*s1*s1) - g * T / s1);
    s2 = s1 + hdt * ps2 / Q;
    //---
    BeW_potential2(pr1, pr2_WW,pr2_BB,pr2_WB, data, U_atom);
    
    for (int i = 0; i < N; i++)
    {
        p1 = data.atoms[i].p;
        f2 = data.atoms[i].f;
        p2 = p1 + hdt * (f2 - (ps2/Q)*p1);
        data.atoms[i].p = p2;
        //s2 = s1 + hdt * ps2 / Q;
    }

    //***new***
    data.s0 = s2;
    data.ps0 = ps2;
    //***new***

    energy(data, pr1, U_atom);
}

void BeW_evolution2(const parameter1 &pr1, const parameter2 &pr2_WW,const parameter2 &pr2_BB,const parameter2 &pr2_WB, Data &data,vector<double> &U_atom)
{
    int N = data.n;
    double dt = pr1.dt;
    double hdt = dt / 2.0;
    //double epsilon = pr1.epsilon;
    double mw1 = 1.0 / pr1.mw;
    double mb1 = 1.0 / pr1.mb;
    //double T = pr1.T;
    //double sigma = pr1.sigma;
    //double steps = pr1.steps;
    //!*****new******
    //double kb = pr1.kb;
    double Q = pr1.Q;
    double g = 3.0 * N;
    double K = 0.0;
    for (int i = 0; i < data.n; i++)
    {
        if(data.atoms[i].name == "W")
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 / pr1.mw);
        else if(data.atoms[i].name == "Be")
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                        + data.atoms[i].p.a10 * data.atoms[i].p.a10
                        + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                        * (0.5 / pr1.mb);
    }
    //double T = data.T;
    double T = pr1.T / pr1.T0;
    double dxi = 2.0 / Q *(K - 3.0/2.0 * N * T);
    //!*****new******
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
        r0 = data.atoms[i].r;
        p0 = data.atoms[i].p;
        f0 = data.atoms[i].f;
        //***new***
        if (data.atoms[i].name == "W")
        {
            p1 = p0 + hdt * f0;
            r2 = r0 + dt * p1 * mw1;
        }
        else if(data.atoms[i].name == "Be")
        {
            p1 = p0 + hdt * f0;
            r2 = r0 + dt * p1 * mb1;
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

    BeW_potential2(pr1, pr2_WW,pr2_BB,pr2_WB, data, U_atom);

    for (int i = 0; i < N; i++)
    {
        p1 = data.atoms[i].p;
        f2 = data.atoms[i].f;
        p2 = p1 + dt * f2 * 0.5;
        data.atoms[i].p = p2;
    }
    energy(data, pr1, U_atom);
}
