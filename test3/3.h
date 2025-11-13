#pragma once

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>


using namespace Eigen;
using namespace std;

// ----------------------------

struct FirstData
{
    int n;                              // number of atoms
    double T;                           // unknown
    double E;                           // energy
    double Box_Ln, Box_Lx,Box_Ly,Box_Lz;// box size
    string name;
    double x, y, z;                     //position
    double vx, vy, vz;                  //velocity
    double dvx, dvy, dvz;               //acceleration
};

struct Data
{
    int n;                              // number of atoms
    double T;                           // unknown
    double E;                           // energy
    double Box_Ln, Box_Lx,Box_Ly,Box_Lz;// box size
    string name;
    double x, y, z;                     //position
    double vx, vy, vz;                  //velocity
    double dvx, dvy, dvz;               //acceleration
};

struct BasicData
{
    string BasicData_filename,Data_filename,Ut_file,Kt_file;
};
