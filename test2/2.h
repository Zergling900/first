#pragma once

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
//#include <eigen3/Eigen/Dense>


//using namespace Eigen;
using namespace std;
// ----------------------------
//int n;        // number of atoms
struct BasicData
{
    double T;     // unknown
    double E;     // energy
    double Box_Ln, Box_Lx,Box_Ly,Box_Lz;// box size
    string data_file_name;//
    string AtomName;//
    double _cell_length_a, _cell_length_b, _cell_length_c;// unit cell sizs
    double _cell_angle_alpha, _cell_angle_beta, _cell_angle_gamma; //unit cell angle
    //double axx,axy,axz,ayx,ayy,ayz,aza,azy,azz;    //
    double randommultiplier;//
};
//----------------------------------------
//  _atom_site_type_symbol
//  _atom_site_label
//  _atom_site_symmetry_multiplicity
//  _atom_site_fract_x
//  _atom_site_fract_y
//  _atom_site_fract_z
//  _atom_site_occupancy
//----------------------------------------


//----------------------------------------
//Matrix3 calculate
//----------------------------------------
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
struct Atom 
{
    string type;
    string label;
    string symmetry_multiplicity;
    double x0, y0, z0;
    string occupancy;
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
    string name;
    double x, y, z;
    double vx, vy, vz;
    double dvx, dvy, dvz;
};

// struct OutPutData
// {
//     string name;
//     double x, y, z;
//     double vx, vy, vz;
//     double dvx, dvy, dvz;
// };

Matrix33 CellVector_L0(const BasicData &data);
Matrix33 CellVector_A(const BasicData &data);
Matrix33 CellVector_C(const Matrix33 &Cell_L0, const Matrix33 &Cell_A);
// Vector3d random(const BasicData &data,const Matrix33 &Cell);

void read1(BasicData &data, vector<Atom> &atoms);
void read2(Data &data);
void build( const vector<Atom> &atoms,
            const Data &data1,
            vector<Data> &datas);
void calculate(int &n,
               const BasicData &data,
               const Matrix33 &Cell,
               const Matrix33 &Cell_L0,
               const vector<Data> &datas,
               vector<Data> &output);

void Output(const int &n, const BasicData &data,
                vector<Data> &output,const Matrix33 &Cell_L0,
                const Matrix33 &Cell,const string &filename);