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
    double Box_Ln, Box_Lx,Box_Ly,Box_Lz;// box size
    string AtomName;//
    double Cell_La, Cell_Lb, Cell_Lc;// unit cell sizs
    double cell_angle_alpha, cell_angle_beta, cell_angle_gamma; //unit cell angle
    double axx,axy,axz,ayx,ayy,ayz,aza,azy,azz;    //
    double randommultiplier;//
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
