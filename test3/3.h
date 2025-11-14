#pragma once

#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

//using namespace Eigen;
using namespace std;

// ----------------------------
//matrix
struct Matrix33
{
    double
    a00,a01,a02,
    a10,a11,a12,
    a20,a21,a22;
    Matrix33()
    : a00(0), a01(0), a02(0),
      a10(0), a11(0), a12(0),
      a20(0), a21(0), a22(0)
{}
    Matrix33(double m00, double m01, double m02,
             double m10, double m11, double m12,
             double m20, double m21, double m22)
        : a00(m00), a01(m01), a02(m02),
          a10(m10), a11(m11), a12(m12),
          a20(m20), a21(m21), a22(m22) {}
};

struct Matrix31
{
    double
    a00,a10,a20;
    Matrix31()
    : a00(0.0), a10(0.0), a20(0.0) {}

    Matrix31(double m00, double m01, double m02)
        : a00(m00), a10(m01), a20(m02) {}
};
//-----------------------------
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
    Matrix33 Box;// box size
    string name;
    double x, y, z;                     //position
    double vx, vy, vz;                  //velocity
    double dvx, dvy, dvz;               //acceleration
};

struct FileName
{
    string BasicData_filename,Data_filename,Ut_file,Kt_file;
    string parameter;
};

struct parameter1
{
    double dt;
    int steps,steps_space;
};

void readF(FileName &filename);