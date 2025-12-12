#include <iostream>
#include <cmath>
#include <vector>

#include "3.h"

void energy(Data &data, const parameter1 &pr1, vector<double> &U_atom)
{
    double m1 = 1.0 / pr1.m;
    data.f_all = sqrt(data.F_all.a00 * data.F_all.a00 + data.F_all.a10 * data.F_all.a10 + data.F_all.a20 * data.F_all.a20);
    data.U_all = 0.0;
    data.K_all = 0.0;
    data.T = 0.0;
    for (int i = 0; i < data.n; i++)
    {
        data.U_all += U_atom[i];
        data.K_all += data.atoms[i].p.a00 * data.atoms[i].p.a00
                    + data.atoms[i].p.a10 * data.atoms[i].p.a10
                    + data.atoms[i].p.a20 * data.atoms[i].p.a20;
    }
    data.K_all = 0.5 * m1 * data.K_all;
    data.T = (2.0 / 3.0) * (data.K_all / (pr1.kb * data.n));
    data.E = data.U_all + data.K_all;
}