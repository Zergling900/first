#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"

using namespace std;

// ----------------------------
// random (first molecular)
std::mt19937 global_eng(std::random_device{}());

void RandomV0(Data &data, const parameter1 &p1)
{
    double sigama0 = sqrt(p1.kb * p1.T * p1.m); // maxwell-boltzmann
    if (sigama0 == 0.0)
    {
        cerr << "[RandomV0] sigma is zero (kb=" << p1.kb << ", T=" << p1.T << ", m=" << p1.m
             << "), velocities will stay 0. Check parameter.p1.\n";
    }
    normal_distribution<double> gauss(0.0, sigama0);//gauss(average, sigma)

    for (int i = 0; i < data.n; i++)
    {
        data.atoms[i].p.a00 = gauss(global_eng);
        data.atoms[i].p.a10 = gauss(global_eng);
        data.atoms[i].p.a20 = gauss(global_eng);
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
}
// ----------------------------

