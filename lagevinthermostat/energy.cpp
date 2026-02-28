#include <iostream>
#include <cmath>
#include <vector>

#include "3.h"

void energy(Data &data, const parameter1 &pr1, vector<double> &U_atom)
{
    data.f_all = sqrt(data.F_all.a00 * data.F_all.a00 + data.F_all.a10 * data.F_all.a10 + data.F_all.a20 * data.F_all.a20);
    data.U_all = 0.0;
    data.K_all = 0.0;
    data.T = 0.0;
    data.P_all = Matrix31(0.0, 0.0, 0.0);
    double s0 = data.s0;
    double ps0 = data.ps0;
    double s02 = s0*s0;
    double ps02 = ps0*ps0;
    double Q = pr1.Q;
    double kb = pr1.kb;
//s
    for (int i = 0; i < data.n; i++)
    {
        data.U_all += U_atom[i];
        data.P_all = data.P_all + data.atoms[i].p ; //--- IGNORE ---
        if(data.atoms[i].name == "W")
            data.K_all += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (1.0 / (2.0 * s02 *pr1.mw));
        else if(data.atoms[i].name == "Be")
            data.K_all += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                        + data.atoms[i].p.a10 * data.atoms[i].p.a10
                        + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                        * (1.0 / (2.0 * s02 *pr1.mb));
    }
    data.T = 2.0 * (data.K_all / (pr1.g * kb));
    data.E = data.U_all + data.K_all;
    data.H = s0*(data.K_all + data.U_all + ps02/(2.0 * Q) + pr1.g * kb * pr1.T * log(s0) - pr1.H0);
}
