#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>//time


#include "3.h"
#include "void.h"



// ----------------------------
// time
// ----------------------------
void BeW_evolution1(const parameter1 &pr1, const parameter2 &pr2_WW,const parameter2 &pr2_BB,const parameter2 &pr2_WB, Data &data,Cell_List &cl,vector<double> &U_atom)
{
auto t0 = std::chrono::high_resolution_clock::now();
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
        const unsigned char type_i = data.atoms[i].atom_type;
        if (type_i == ATOM_TYPE_W)
        {
            r2 = r0 + (dt * p1 * (1.0/ (mw * s1)));
        }
        else if (type_i == ATOM_TYPE_BE)
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
        if (data.atoms[i].atom_type == ATOM_TYPE_W)
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 /(pr1.mw * s1 * s1));
        else if (data.atoms[i].atom_type == ATOM_TYPE_BE)
            K += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 / (pr1.mb * s1 * s1));
    }
    ps3 = ps2 - dt *(-K + g *kb * T * log(s1) - H0 + g *kb * T);
    //***d***
    //---
//auto t1 = std::chrono::high_resolution_clock::now();
    // lcl1(data, cl, pr1, pr2_WW, U_atom);
    // lcl2(data, cl, pr1, pr2_WW, pr2_BB, pr2_WB, U_atom);
    Force_Current(pr1, pr2_WW, pr2_BB, pr2_WB, data, cl, U_atom);
//auto t2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++)
    {
        p1 = data.atoms[i].p;
        f2 = data.atoms[i].f;
        p2 = p1 + hdt * s1 *(f2);
        data.atoms[i].p = p2;
    }
    u2 = 0.0;
    for (int i = 0; i < N; ++i)
    {
        u2 += U_atom[i];
    }
    ps4 = ps3 - hdt * u2;
    //***e***
    s2 = s1 * (1.0 + hdt * ps4 / (2.0 * Q))*(1.0 + hdt * ps4 / (2.0 * Q));
    ps5 = ps4 /(1.0+ hdt * ps4 / (2.0 * Q));


    //***new***
    data.s0 = s2;
    data.ps0 = ps5;
    // Recompute observables with the updated thermostat variables so H/K/T
    // are consistent with the final state of this time step.
    energy(data, pr1, U_atom);
    //***new***
//auto t3 = std::chrono::high_resolution_clock::now();
/*
    double ms0 =
    std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms1 =
    std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::cout << "time_evo = " << ms0+ms1 << " ms\n";
    double ms2 =
    std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "time_potential = " << ms2 << " ms\n";
*/
}

void BeW_evolution2(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,Cell_List &cl,vector<double> &U_atom)
{
    const int N = data.n;
    const double dt = pr1.dt;
    const double hdt = 0.5 * dt;
    const double mw = pr1.mw;
    const double mb = pr1.mb;

    const double xi = pr1.xi;
    if (xi < 0.0)
    {
        std::cerr << "[BeW_evolution2] xi < 0 is invalid for Langevin: " << xi << "\n";
        return;
    }

    static thread_local std::mt19937 rng(std::random_device{}());

    const double c = std::exp(-xi * dt);
    double one_minus_c2 = 1.0 - c * c;
    if (one_minus_c2 < 0.0)
        one_minus_c2 = 0.0;

    const double kbT = pr1.kb * pr1.T;
    const double sigmaW = (kbT > 0.0) ? std::sqrt(mw * kbT * one_minus_c2) : 0.0;
    const double sigmaBe = (kbT > 0.0) ? std::sqrt(mb * kbT * one_minus_c2) : 0.0;
    std::normal_distribution<double> gaussW(0.0, sigmaW);
    std::normal_distribution<double> gaussBe(0.0, sigmaBe);

    auto wrap_coord = [](double x, double L)
    {
        while (x >= L)
            x -= L;
        while (x < 0.0)
            x += L;
        return x;
    };

    // Langevin thermostat does not use Nose extended variables.
    data.s0 = 1.0;
    data.ps0 = 0.0;

    for (int i = 0; i < N; ++i)
    {
        Atom &a = data.atoms[i];
        const bool isW = (a.atom_type == ATOM_TYPE_W);
        const double mass = isW ? mw : mb;

        // BAOAB-like step in momentum form:
        // half kick -> half drift -> OU thermostat -> half drift
        Matrix31 p = a.p + hdt * a.f;
        Matrix31 r = a.r + (hdt / mass) * p;

        p = c * p;
        if (isW)
        {
            if (sigmaW > 0.0)
            {
                p.a00 += gaussW(rng);
                p.a10 += gaussW(rng);
                p.a20 += gaussW(rng);
            }
        }
        else
        {
            if (sigmaBe > 0.0)
            {
                p.a00 += gaussBe(rng);
                p.a10 += gaussBe(rng);
                p.a20 += gaussBe(rng);
            }
        }

        r = r + (hdt / mass) * p;
        r.a00 = wrap_coord(r.a00, data.Box.a00);
        r.a10 = wrap_coord(r.a10, data.Box.a11);
        r.a20 = wrap_coord(r.a20, data.Box.a22);

        a.r = r;
        a.p = p;
    }

    // lcl1(data, cl, pr1, pr2_WW, U_atom);
    // lcl2(data, cl, pr1, pr2_WW, pr2_BB, pr2_WB, U_atom);
    Force_Current(pr1, pr2_WW, pr2_BB, pr2_WB, data, cl, U_atom);

    for (int i = 0; i < N; ++i)
    {
        data.atoms[i].p = data.atoms[i].p + hdt * data.atoms[i].f;
    }

    energy(data, pr1, U_atom);
}
