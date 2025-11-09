#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <eigen3/Eigen/Dense>

#include "2.h"

using namespace Eigen;
using namespace std;

// ----------------------------
Matrix3d CellVector(const BasicData &data)
{
    const double a = data.Cell_La;
    const double b = data.Cell_Lb;
    const double c = data.Cell_Lc;

    //
    const double alpha = data.cell_angle_alpha * M_PI / 180.0;
    const double beta = data.cell_angle_beta * M_PI / 180.0;
    const double gamma = data.cell_angle_gamma * M_PI / 180.0;

    const double ca = std::cos(alpha);
    const double cb = std::cos(beta);
    const double cg = std::cos(gamma);
    const double sg = std::sin(gamma);

    Matrix3d Cell; // 3*3

    // 1 list a = a*(1, 0, 0)
    Cell(0, 0) = data.Cell_La;
    Cell(1, 0) = 0.0;
    Cell(2, 0) = 0.0;

    // 2 list b = b*(cos(γ), sin(γ), 0 )
    Cell(0, 1) = data.Cell_Lb * cg;
    Cell(1, 1) = data.Cell_Lb * sg;
    Cell(2, 1) = 0.0;

    // 3 list c = c* (cos(β), (cos(α) - cos(β)*cos(γ)) / sin(γ),
    //                1 + 2*cos(α)*cos(β)*cos(γ)- cos(α)*cos(α)
    //                - cos(β)*cos(β)- cos(γ)*cos(γ)) / sin(γ))
    Cell(0, 2) = data.Cell_Lc * cb;
    Cell(1, 2) = data.Cell_Lc * (ca - cb * cg) / sg;
    Cell(2, 2) = data.Cell_Lc * std::sqrt(1.0 + 2.0 * ca * cb * cg - ca * ca - cb * cb - cg * cg) / sg;

    return Cell;
}
//-----------------------------
Vector3d Step(const Matrix3d &Cell)
{
    return Vector3d(
        Cell.col(0).norm(),
        Cell.col(1).norm(),
        Cell.col(2).norm());
}

// random (first molecular)
// std::mt19937 global_eng(std::random_device{}());
// Vector3d random(const BasicData &data, const Matrix3d &Cell)
// {
//     double f = data.Cell_La;

//     if (f < data.Cell_Lb)
//         f = data.Cell_Lb;

//     if (f < data.Cell_Lc)
//         f = data.Cell_Lc;
//     std::uniform_real_distribution<> distr(0, 0.0 * f);
//     Vector3d duvw(distr(global_eng), distr(global_eng), distr(global_eng));
//     return Cell * duvw;
// }
//