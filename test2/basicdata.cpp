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
Matrix3d Cell0(const BasicData &data)
{    
    Matrix3d Cell0;
    Cell0 << data.Box_Lx, 0, 0,
             0, data.Box_Ly, 0,
             0, 0, data.Box_Lz;

    Matrix3d R;
        R = AngleAxisd(data.thetayz, Vector3d::UnitX())   //  x
        * AngleAxisd(data.thetaxz, Vector3d::UnitY())   //  y
        * AngleAxisd(data.thetaxy, Vector3d::UnitZ());  //  z
    return R * Cell0;
}

//-----------------------------
Vector3d Step(const Matrix3d &Cell)
{
    return Vector3d(
        Cell.col(0).norm(),
        Cell.col(1).norm(),
        Cell.col(2).norm()
    );
}


// random (first molecular)
std::mt19937 global_eng(std::random_device{}());
void random(const BasicData &data, FirstMolecularData &MolecularData, double scale, const Matrix3d &Cell)
{
double f = data.Cell_Lz;

if (f < data.Cell_Ly)
    f = data.Cell_Ly;

if (f < data.Cell_Lx)
    f = data.Cell_Lx;
    std::uniform_real_distribution<> distr(-scale, scale); 
    Vector3d duvw(distr(global_eng), distr(global_eng), distr(global_eng));
    Vector3d r0 = Cell * duvw;

    MolecularData.x0 = r0(0);
    MolecularData.y0 = r0(1);
    MolecularData.z0 = r0(2);
}
//

