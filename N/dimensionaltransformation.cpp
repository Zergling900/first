#include <iostream>
#include <cmath>
#include <vector>

#include "3.h"

void Dimensionalless(Data &data, parameter1 &pr1, parameter2 &pr2_WW, parameter2 &pr2_BB, parameter2 &pr2_WB)
{
    double tau = pr1.tau;
    //V = a \cdot (b \times c)
    pr1.length = cbrt((data.Box.a00 * (data.Box.a11 * data.Box.a22 - data.Box.a12 * data.Box.a21)
                  + data.Box.a01 * (data.Box.a12 * data.Box.a20 - data.Box.a10 * data.Box.a22)
                  + data.Box.a02 * (data.Box.a10 * data.Box.a21 - data.Box.a11 * data.Box.a20))/data.n);
    double length = pr1.length;
    //box length
    data.Box = data.Box * (1.0/length);

    //mass
    double m = pr1.mb0;
    pr1.E0 = (m*length*length)/(tau*tau)*103.642695;//eV
    pr1.F0 = (m * length / (tau * tau));
    pr1.P0 = (m * length / tau);
    pr1.T0 = pr1.E0 / pr1.kb;

    pr1.mw = pr1.mw0 / m;
    pr1.mb = pr1.mb0 / m; 
    pr1.T = pr1.T / pr1.T0;
    pr1.Q = pr1.Q/(pr1.E0 * pr1.tau * pr1.tau);
    //atom
    for (int i = 0; i < data.n; i++)
    {
        data.atoms[i].r = data.atoms[i].r *(1.0/length);
    }

    //time
    pr1.dt = pr1.dt / tau;
    //pr1.endtime = pr1.endtime / pr1.tau;

    pr2_WW.D0 = pr2_WW.D0 / pr1.E0;
    pr2_BB.D0 = pr2_BB.D0 / pr1.E0;
    pr2_WB.D0 = pr2_WB.D0 / pr1.E0;

    pr2_WW.r0 = pr2_WW.r0 / length;
    pr2_BB.r0 = pr2_BB.r0 / length;
    pr2_WB.r0 = pr2_WB.r0 / length;

    pr2_WW.beta = pr2_WW.beta * length;
    pr2_BB.beta = pr2_BB.beta * length;
    pr2_WB.beta = pr2_WB.beta * length;

    pr2_WW.R = pr2_WW.R / length;
    pr2_BB.R = pr2_BB.R / length;
    pr2_WB.R = pr2_WB.R / length;

    pr2_WW.D = pr2_WW.D / length;
    pr2_BB.D = pr2_BB.D / length;
    pr2_WB.D = pr2_WB.D / length;

    //！test （this right？）
    pr2_WW.mu = pr2_WW.mu * length;
    pr2_BB.mu = pr2_BB.mu * length;
    pr2_WB.mu = pr2_WB.mu * length;

}