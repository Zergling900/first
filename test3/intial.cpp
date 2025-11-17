#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"

using namespace std;

// ----------------------------
// random (first molecular)
std::mt19937 global_eng(std::random_device{}());

void random( FirstData&Data0)
{
    std::uniform_real_distribution<> distr();
    Data0.vx = distr(global_eng);
    Data0.vy = distr(global_eng);
    Data0.vz = distr(global_eng);

}