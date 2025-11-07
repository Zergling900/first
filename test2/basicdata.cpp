#pragma once

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <eigen3/Eigen/Dense>


using namespace Eigen;
using namespace std;
// ----------------------------
struct BasicData
{
    int n;        // number of atoms
    double T;     // unknown
    double E;     // energy
    double Box_Lx;// box size
    double Box_Ly;// box size
    double Box_Lz;// box size
    double a0;    //
    double theta; //
    double phi;   //
};

struct FirstMolecularData
{
    double x0, y0, z0;
};

//-----------------------------------------------
// struct AtomicData
// {
//     double x1, y1, z1;
//     double xv1, yv1, zv1;
//     double xdv1, ydv1, zdv1;
// };
//std::vector<AtomicData> BCC;
// std::vector<AtomicData> FCC;
// std::vector<AtomicData> Diamond;
// -----------------------------------------------
struct Data
{
    double x, y, z;
    double xv, yv, zv;
    double xdv, ydv, zdv;
};
using BCCData = Data;
using FCCData = Data;
using DiamondData = Data;
std::vector<BCCData> BCC;
std::vector<FCCData> FCC;
std::vector<DiamondData> Diamond;

Matrix3d A0(const double &a0)
{
    Matrix3d A;
    A << data.a0, 0, 0,
        0, data.a0, 0,
        0, 0, data.a0;
    return A;
};

Matrix3d BOX(const BasicData &data) 
{
    Matrix3d Box_L;
    Box_L << data.Box_Lx, 0, 0,
            0, data.Box_Ly, 0,
            0, 0, data.Box_Lz;
    return Box_L;
};
// random (first molecular)
std::mt19937 global_eng(std::random_device{}());
void random(FirstMolecularData &MolecularData, double a0)
{
    std::uniform_real_distribution<> distr(0.0, 0.1 * a0);
    MolecularData.x0 = distr(global_eng);
    MolecularData.y0 = distr(global_eng);
    MolecularData.z0 = distr(global_eng);
};
//
// Rotation matrix (theta, phi)
Matrix3d Rotation(double theta, double phi)
{
    double ct = std::cos(theta);
    double st = std::sin(theta);
    double cp = std::cos(phi);
    double sp = std::sin(phi);

    Matrix3d Rz;
    Rz << cp, -sp, 0,
        sp, cp, 0,
        0, 0, 1;

    Matrix3d Ry;
    Ry << ct, 0, st,
        0, 1, 0,
        -st, 0, ct;

    Matrix3d R = Ry * Rz;
    return R;
};

double Step(const Vector3d &a)
{
    double ax = std::abs(a(0));
    double ay = std::abs(a(1));
    double az = std::abs(a(2));
    return std::max({ax, ay, az});
}
