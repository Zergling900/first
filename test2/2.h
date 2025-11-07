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
    double Box_Lx,Box_Ly,Box_Lz;// box size
    double Cell_Lx, Cell_Ly, Cell_Lz;// unit cell sizs
    double axx,axy,axz,ayx,ayy,ayz,aza,azy,azz;    //

    double thetaxy,thetaxz,thetayz; //
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
