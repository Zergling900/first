#pragma once

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <eigen3/Eigen/Dense>


void Box(const Matrix3d &A_rot, const BasicData &data, int &nx_est, int &ny_est, int &nz_est)
{
    // col(0) = (a0x, a0y, a0z)
    Vector3d a1 = A_rot.col(0);
    Vector3d a2 = A_rot.col(1);
    Vector3d a3 = A_rot.col(2);

    double dx = Step(a1);
    double dy = Step(a2);
    double dz = Step(a3);

    if (dx < 1e-6)
        dx = 1e-6;
    if (dy < 1e-6)
        dy = 1e-6;
    if (dz < 1e-6)
        dz = 1e-6;

    nx_est = static_cast<int>(std::ceil(data.Box_L / dx)) + 1;
    ny_est = static_cast<int>(std::ceil(data.Box_L / dy)) + 1;
    nz_est = static_cast<int>(std::ceil(data.Box_L / dz)) + 1;
}