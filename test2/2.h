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

Matrix3d CellVector_L0(const BasicData &data);
Matrix3d CellVector_A(const BasicData &data);
Matrix3d CellVector_C(const Matrix3d &Cell_L0, const Matrix3d &Cell_A);
// Vector3d random(const BasicData &data,const Matrix3d &Cell);

void read1(BasicData &data, vector<Atom> &atoms);
void read2(Data &data);
void build( const vector<Atom> &atoms,
            const Data &data1,
            vector<Data> &datas);
void calculate(int &n,
               const BasicData &data,
               const Matrix3d &Cell,
               const vector<Data> &datas,
               vector<Data> &output);

void Output(const int &n, const BasicData &data,
                vector<Data> &output,const Matrix3d &Cell,
                const string &filename);