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
    double s2 = data.s0 *data.s0;
    for (int i = 0; i < data.n; i++)
    {
        data.U_all += U_atom[i];
        if(data.atoms[i].name == "W")
            data.K_all += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                    * (0.5 / (pr1.mw *s2));
        else if(data.atoms[i].name == "Be")
            data.K_all += (data.atoms[i].p.a00 * data.atoms[i].p.a00
                        + data.atoms[i].p.a10 * data.atoms[i].p.a10
                        + data.atoms[i].p.a20 * data.atoms[i].p.a20)
                        * (0.5 / (pr1.mb *s2));
    }
    data.T = 2.0 * (data.K_all / pr1.g);
    data.E = data.U_all + data.K_all;
    data.H = data.E + data.ps0 * data.ps0 / (2.0 * pr1.Q) + pr1.g * pr1.T * log(data.s0);
}