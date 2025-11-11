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
Matrix3d CellVector_L0(const BasicData &data)
{
    const double a = data._cell_length_a;
    const double b = data._cell_length_b;
    const double c = data._cell_length_c;

    Matrix3d Cell_L0; // 3*3
    Cell_L0(0, 0) = data._cell_length_a;
    Cell_L0(1, 0) = 0.0;
    Cell_L0(2, 0) = 0.0;

    // 2 list b = b*(cos(γ), sin(γ), 0 )
    Cell_L0(0, 1) = 0.0;
    Cell_L0(1, 1) = data._cell_length_b;
    Cell_L0(2, 1) = 0.0;

    // 3 list c = c* (cos(β), (cos(α) - cos(β)*cos(γ)) / sin(γ),
    //                1 + 2*cos(α)*cos(β)*cos(γ)- cos(α)*cos(α)
    //                - cos(β)*cos(β)- cos(γ)*cos(γ)) / sin(γ))
    Cell_L0(0, 2) = 0.0;
    Cell_L0(1, 2) = 0.0;
    Cell_L0(2, 2) = data._cell_length_c;

    return Cell_L0;
}

Matrix3d CellVector_A(const BasicData &data)
{
    const double alpha = data._cell_angle_alpha * M_PI / 180.0;
    const double beta = data._cell_angle_beta * M_PI / 180.0;
    const double gamma = data._cell_angle_gamma * M_PI / 180.0;

    const double ca = std::cos(alpha);
    const double cb = std::cos(beta);
    const double cg = std::cos(gamma);
    const double sg = std::sin(gamma);

    Matrix3d Cell_A; // 3*3
    // 1 list a = a*(1, 0, 0)
    Cell_A(0, 0) = 1;
    Cell_A(1, 0) = 0.0;
    Cell_A(2, 0) = 0.0;

    // 2 list b = b*(cos(γ), sin(γ), 0 )
    Cell_A(0, 1) = cg;
    Cell_A(1, 1) = sg;
    Cell_A(2, 1) = 0.0;

    // 3 list c = c* (cos(β), (cos(α) - cos(β)*cos(γ)) / sin(γ),
    //                1 + 2*cos(α)*cos(β)*cos(γ)- cos(α)*cos(α)
    //                - cos(β)*cos(β)- cos(γ)*cos(γ)) / sin(γ))
    Cell_A(0, 2) = cb;
    Cell_A(1, 2) = (ca - cb * cg) / sg;
    Cell_A(2, 2) = std::sqrt(1.0 + 2.0 * ca * cb * cg - ca * ca - cb * cb - cg * cg) / sg;

    return Cell_A;
}

Matrix3d CellVector_C(const Matrix3d &Cell_L0, const Matrix3d &Cell_A)
{
    Matrix3d Cell;
    Cell = Cell_A * Cell_L0;

    return Cell;
}
//-----------------------------
// random (first molecular)
// std::mt19937 global_eng(std::random_device{}());
// Vector3d random(const BasicData &data, const Matrix3d &Cell)
// {
//     double f = data._cell_length_a;

//     if (f < data._cell_length_b)
//         f = data._cell_length_b;

//     if (f < data._cell_length_c)
//         f = data._cell_length_c;
//     std::uniform_real_distribution<> distr(0, 0.0 * f);
//     Vector3d duvw(distr(global_eng), distr(global_eng), distr(global_eng));
//     return Cell * duvw;
// }
//