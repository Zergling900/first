#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

#include "3.h"
#include "void.h"
using namespace std;

struct parameter2
{
    double dt, T, kb, m;
    int steps, steps_space;
    double D0, r0, beta, S;
    double gamma, c, d, h;
    double R, D;
    double mu, rf, bf;
};

void read2(const FileName &filename, parameter2 &pr2)
{
    ifstream fin(filename.parameter2_filename);
    if (!fin)
    {
        std::cerr << "can't open file: " << filename.parameter2_filename << std::endl;
        return;
    }

    string key, line;
    double value;
    bool reading_abop = false;

    while (getline(fin, line))
    {
        // 去掉注释
        size_t pos = line.find("//");
        if (pos != string::npos)
            line = line.substr(0, pos);

        stringstream ss(line);

        // ----------------------------------------------------------------

        if (key == "ABOP")
        {
            ss >> key; // Abb / Aww / Abw
            cout << "Using ABOP parameter set: " << key << endl;
            reading_abop = true;
        }
        else if (reading_abop)
        {
            if (key == "D0")
                ss >> pr2.D0;
            else if (key == "r0")
                ss >> pr2.r0;
            else if (key == "beta")
                ss >> pr2.beta;
            else if (key == "S")
                ss >> pr2.S;
            else if (key == "gamma")
                ss >> pr2.gamma;
            else if (key == "c")
                ss >> pr2.c;
            else if (key == "d")
                ss >> pr2.d;
            else if (key == "h")
                ss >> pr2.h;
            else if (key == "R")
                ss >> pr2.R;
            else if (key == "D")
                ss >> pr2.D;
            else if (key == "mu2")
                ss >> pr2.mu;
            else if (key == "rf")
                ss >> pr2.rf;
            else if (key == "bf")
                ss >> pr2.bf;
        }
    }

    fin.close();
};

void BeWpotential(const parameter1 &pr1, const parameter2 &pr2, Data &data, vector<double> &U_atom)
{

    int i, j, k;
    // int N = data.n;
    // data.U_all = 0.0;
    // data.K_all = 0.0;
    data.F_all = Matrix31(0.0, 0.0, 0.0);
    data.f_all = 0.0;

    double D0 = pr2.D0;
    double R = pr2.R;
    double D = pr2.D;
    double mu = pr2.mu;
    double rf = pr2.rf;
    double bf = pr2.bf;
    double r0 = pr2.r0;
    double beta = pr2.beta;
    double S = pr2.S;
    double gamma = pr2.gamma;
    double c = pr2.c;
    double d = pr2.d;
    double h = pr2.h;

    double m = pr1.m;

    U_atom.assign(data.n, 0.0); // initialize

    Matrix31 F = Matrix31(0.0, 0.0, 0.0);
    //--------------------------------------------------------
    double Lx = data.Box.a00;
    double Ly = data.Box.a11;
    double Lz = data.Box.a22;
    double Lxh = data.Box.a00 / 2.0;
    double Lyh = data.Box.a11 / 2.0;
    double Lzh = data.Box.a22 / 2.0;

    // double rc2 = r0 * r0;

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
            Matrix31 drij = aj.r - ai.r;                                                     // distance
            double rij2 = (drij.a00 * drij.a00 + drij.a10 * drij.a10 + drij.a20 * drij.a20); // distance
            double rij = sqrt(rij2);                                                         // distance

            double VR = D0 / (S - 1) * exp(-beta * sqrt(2 * S) * (rij - r0));
            double VA = S * VR;

            double fc1 = 1.0;                                         // r <= R-D
            double fc2 = 0.5 - 0.5 * sin(0.5 * M_PI * (rij - R) / D); // |R-dr| <= D
            double fc3 = 0.0;                                         // r >= R+D

            double Xij = 0.0;
            double Xji = 0.0;

            for (k = 0; k < data.n; k++)
            {
                if (k == i || k == j)
                    continue;

                const Atom &ak = data.atoms[k];
                Matrix31 drik = ak.r - ai.r; // distance
                Matrix31 drjk = ak.r - aj.r;

                double rik2 = (drik.a00 * drik.a00 + drik.a10 * drik.a10 + drik.a20 * drik.a20); // distance
                double rjk2 = (drjk.a00 * drjk.a00 + drjk.a10 * drjk.a10 + drjk.a20 * drjk.a20); // distance

                double rik = sqrt(rik2);                                   // distance
                double rjk = sqrt(rjk2);                                   // distance
                double fck1 = 1.0;                                         // r <= R-D
                double fck2 = 0.5 - 0.5 * sin(0.5 * M_PI * (rik - R) / D); // |R-dr| <= D
                double fck22 = 0.5 - 0.5 * sin(0.5 * M_PI * (rjk - R) / D);
                double fck3 = 0.0;                                         // r >= R+D
                //------------------------------------------------------------
                // bji
                //------------------------------------------------------------
                double cost_ijk = (rij2 + rik2 - rjk2) / (2.0 * rij * rik); // theta_ijk angle is i?
                
                if (cost_ijk > 1.0)
                    cost_ijk = 1.0;
                if (cost_ijk < -1.0)
                    cost_ijk = -1.0;

                double g_ijk = gamma * (1 + c * c * (1.0 / (d * d) - 1 / (d * d + (h + cost_ijk) * (h + cost_ijk))));

                if (rik <= R - D)
                {
                    Xij += fck1 * g_ijk * exp(mu * (rij - rik));
                }
                else if (rik >= R + D)
                {
                    Xij += 0.0;
                }
                else
                {
                    Xij += fck2 * g_ijk * exp(mu * (rij - rik));
                }

                //------------------------------------------------------------
                // bji
                //------------------------------------------------------------
                double cost_jik = (rij2 + rjk2 - rik2) / (2.0 * rij * rjk); // theta_ijk angle is j?
                
                if (cost_jik > 1.0)
                    cost_jik = 1.0;
                if (cost_jik < -1.0)
                    cost_jik = -1.0;

                double g_jik = gamma * (1 + c * c * (1.0 / (d * d) - 1 / (d * d + (h + cost_jik) * (h + cost_jik))));

                if (rjk <= R - D)
                {
                    Xji += fck1 * g_jik * exp(mu * (rij - rjk));
                }
                else if (rjk >= R + D)
                {
                    Xji += 0.0;
                }
                else
                {
                    Xji += fck22 * g_jik * exp(mu * (rij - rjk));
                }
            }
            double bij = 1.0 / sqrt(1 + Xij);
            double bji = 1.0 / sqrt(1 + Xji);

        }

        
    }
}