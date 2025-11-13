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


