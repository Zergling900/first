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
               const Matrix3d &Cell_L0,
               const vector<Data> &datas,
               vector<Data> &output)
{
    // --------------------------------------------
    output.clear();

    for (int ix = 0; ix < data.Box_Ln + 1; ++ix)
    {
        for (int iy = 0; iy < data.Box_Ln + 1; ++iy)
        {
            for (int iz = 0; iz < data.Box_Ln + 1; ++iz)
            {
                for (const auto &d : datas)
                {
                    Data aaa = d;

                    Vector3d frac(d.x, d.y, d.z);
                    Vector3d step(ix, iy, iz);
                    Vector3d pos = frac + step;
                    //Vector3d r = pos;
                    Vector3d r = Cell_L0 * pos;
                    aaa.x = r.x();
                    aaa.y = r.y();
                    aaa.z = r.z();

                    output.push_back(aaa);
                }
            }
        }
    }
    n = static_cast<int>(output.size());
}

void Output(const int &n, const BasicData &data,
                vector<Data> &output,const Matrix3d &Cell_L0,
                const Matrix3d &Cell,const string &filename)
{
    FILE *fp = fopen(filename.c_str(), "w");
    Matrix3d Box = Cell_L0 * (data.Box_Ln + 1.0);
    //Matrix3d Box = Cell * (data.Box_Ln + 1.0);
    fprintf(fp, "%d\n", n);
    fprintf(fp, "   time=   %f (fs)  Energy=  %f (eV)\n", data.T, data.E);
    fprintf(fp, "BOX %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
           Box(0, 0), Box(1, 0), Box(2, 0), Box(0, 1), Box(1, 1), Box(2, 1), Box(0, 2), Box(1, 2), Box(2, 2));
    for(const auto &d : output)
    {
        fprintf(fp, " %4s  %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            d.name.c_str(),d.x,d.y,d.z,d.vx,d.vy,d.vz,d.dvx,d.dvy,d.dvz);
    }
    fclose(fp);
}