#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
//#include <eigen3/Eigen/Dense>

#include "2.h"

//using namespace Eigen;
using namespace std;

//----------------------------------------
// make first unit data
void build(const vector<Atom> &atoms,
           const Data &data1,
           vector<Data> &datas)
{
    datas.clear();
    datas.reserve(atoms.size());

    for (const auto &atom : atoms)
    {
        Data d = data1;
        d.name = atom.type;
        d.x = atom.x0;
        d.y = atom.y0;
        d.z = atom.z0;
        datas.push_back(d);
    }
}
void calculate(int &n,
               const BasicData &data,
               const Matrix33 &Cell,
               const Matrix33 &Cell_L0,
               const vector<Data> &datas,
               vector<Data> &output)
{
    // --------------------------------------------
    output.clear();

    for (int ix = 0; ix < data.Box_Ln + 1 ; ++ix)
    {
        for (int iy = 0; iy < data.Box_Ln+ 1 ; ++iy)
        {
            for (int iz = 0; iz < data.Box_Ln + 1 ; ++iz)
            {
                for (const auto &d : datas)
                {
                    Data aaa = d;

                    Matrix31 frac(d.x, d.y, d.z);
                    Matrix31 step(ix, iy, iz);
                    Matrix31 pos = frac + step;
                    //Matrix31 r = pos;
                    Matrix31 r = Cell * pos;
                    //Matrix31 r = Cell_L0 * pos;
                    aaa.x = r.a00;
                    aaa.y = r.a10;
                    aaa.z = r.a20;

                    output.push_back(aaa);
                }
            }
        }
    }
    n = static_cast<int>(output.size());
}

void Output(const int &n, const BasicData &data,
                vector<Data> &output,const Matrix33 &Cell_L0,
                const Matrix33 &Cell,const string &filename)
{
    FILE *fp = fopen(filename.c_str(), "w");
    //Matrix33 Box = Cell_L0 * (data.Box_Ln + 1.0);
    Matrix33 Box = Cell * (data.Box_Ln + 1);
    fprintf(fp, "%d\n", n);
    fprintf(fp, "   time=   %f (fs)  Energy=  %f (eV)\n", data.T, data.E);
    fprintf(fp, "BOX %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
         Box.a00, Box.a10, Box.a20, Box.a01, Box.a11, Box.a21, Box.a02, Box.a12, Box.a22);
    for(const auto &d : output)
    {
        fprintf(fp, " %4s  %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            d.name.c_str(),d.x,d.y,d.z,d.vx,d.vy,d.vz,d.dvx,d.dvy,d.dvz);
    }
    fclose(fp);
}