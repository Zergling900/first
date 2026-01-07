#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"



// ----------------------------
// random (first molecular)
std::mt19937 global_eng(std::random_device{}());

void RandomV0(Data &data, parameter1 &pr1)
{
    double sigamaW = std::sqrt(pr1.kb * pr1.T * pr1.mw); // maxwell-boltzmann(W)
    double sigamaBe = std::sqrt(pr1.kb * pr1.T * pr1.mb); // maxwell-boltzmann(Be)
    if (sigamaW == 0.0 || sigamaBe == 0.0)
    {
        cerr << "[RandomV0] sigma is zero (kb=" << pr1.kb << ", T=" << pr1.T << ", mW=" << pr1.mw << ", mBe=" << pr1.mb
             << "), velocities will stay 0. Check parameter.pr1.\n";
    }
    normal_distribution<double> gaussW(0.0, sigamaW);//gauss(average, sigma)
    normal_distribution<double> gaussBe(0.0, sigamaBe);//gauss(average, sigma)

    for (int i = 0; i < data.n; i++)
    {
        if (data.atoms[i].name == "W")
        {
            data.atoms[i].p.a00 = gaussW(global_eng);
            data.atoms[i].p.a10 = gaussW(global_eng);
            data.atoms[i].p.a20 = gaussW(global_eng);
        }
        else if (data.atoms[i].name == "Be")
        {
            data.atoms[i].p.a00 = gaussBe(global_eng);
            data.atoms[i].p.a10 = gaussBe(global_eng);
            data.atoms[i].p.a20 = gaussBe(global_eng);
        }
       
    }

    //p_cm = 0
    Matrix31 p_cm(0.0, 0.0, 0.0);

    for (int i = 0; i < data.n; i++)
    {
        p_cm = data.atoms[i].p + p_cm;
    }
    //
    p_cm.a00 /= data.n;
    p_cm.a10 /= data.n;
    p_cm.a20 /= data.n;

    //
    for (int i = 0; i < data.n; i++)
    {
        data.atoms[i].p = data.atoms[i].p - p_cm;
    }
    pr1.g = 3 * data.n - 3;
    data.s0 = pr1.s0;
    data.ps0 = pr1.ps0;
}
// ----------------------------

