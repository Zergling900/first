#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"

using namespace std;

// ----------------------------
// random (first molecular)
std::mt19937 global_eng(std::random_device{}());

void RandomV0(Data &Data0, const parameter1 &p1)
{
    double sigama0 = sqrt(p1.kb * Data0.T / p1.m); // maxwell-boltzmann
    std::normal_distribution<double> gauss(0.0, sigama0);//gauss(average, sigma)

    for (int i = 0; i < Data0.n; i++)
    {
        Data0.atoms[i].v.a00 = gauss(global_eng);
        Data0.atoms[i].v.a10 = gauss(global_eng);
        Data0.atoms[i].v.a20 = gauss(global_eng);
    }
    //V_cm = 0
    Matrix31 v_cm(0, 0, 0);

    for (int i = 0; i < Data0.n; i++)
    {
        v_cm.a00 += Data0.atoms[i].v.a00;
        v_cm.a10 += Data0.atoms[i].v.a10;
        v_cm.a20 += Data0.atoms[i].v.a20;
    }

    v_cm.a00 /= Data0.n;
    v_cm.a10 /= Data0.n;
    v_cm.a20 /= Data0.n;

    //
    for (int i = 0; i < Data0.n; i++)
    {
        Data0.atoms[i].v.a00 -= v_cm.a00;
        Data0.atoms[i].v.a10 -= v_cm.a10;
        Data0.atoms[i].v.a20 -= v_cm.a20;
    }
}
// ----------------------------



