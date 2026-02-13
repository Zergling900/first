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

//---------------------------------------------------------------------------------------------------------
//Matrix calculate
//33 * 33
Matrix33 operator*(const Matrix33 &A,const Matrix33 &B);
//33 +33
Matrix33 operator+(const Matrix33 &A,const Matrix33 &B);
//33 -33
Matrix33 operator-(const Matrix33 &A,const Matrix33 &B);
//33 *31
Matrix31 operator*(const Matrix33 &A,const Matrix31 &B);
//31 +31
Matrix31 operator+(const Matrix31 &A,const Matrix31 &B);
//31 - 31
Matrix31 operator-(const Matrix31 &A,const Matrix31 &B);
//k*31
Matrix31 operator*(double k, const Matrix31 &v);
//31*k
Matrix31 operator*(const Matrix31 &v, double k);
//k*33
Matrix33 operator*(double k, const Matrix33 &M);
//33*k
Matrix33 operator*(const Matrix33 &M, double k);
//---------------------------------------------------------------------------------------------------------

// struct AtomPotential
// {
    
// };
struct Atom
{
    string name;
    Matrix31 r;                     //position                    //velocity
    Matrix31 p;                     //momentum
    Matrix31 f;                     //force
};

struct Cell_List
{
    int Mx,My,Mz;
    double Wx,Wy,Wz;
    int cell_num;
    std::vector<int> Cell;
    std::vector<int> cell_offset;
    std::vector<int> atom_indices;
};

struct Data
{
    int n;                              // number of atoms
    double t,T,E,H,s0,ps0;                           // unknown
    double U_all,K_all,f_all;                           // energy
    Matrix31 F_all,P_all;
    Matrix33 Box;
    std::vector<Atom> atoms;                       //momentum
};

struct FileName
{
    string BasicData_filename,Data_filename,Et_file;
    string parameter_filename,parameter2_filename,parameter3_filename;
};

struct parameter1
{
    double dt, epsilon, kb,T, TT,T_end, sigma, mw, mb, endtime;
    double output_Data_time, output_Et_time, T_time, dT;
    double s0,ps0, xi, Q;
    double E0,H0;
    int g;
    int steps,steps_space;
};

struct parameter2
{
    double D0, r0, beta, S;
    double gamma, c, d, h;
    double R, D;
    double mu, rf, bf;
};

struct parameter3
{
    int extra_steps_max = 400000;
    int plateau_blocks = 1;
    int output_rotate_every_temp_up = 10;
    int output_rotate_every_temp_down = 10;
    int output_flush_every_writes = 256;
};
