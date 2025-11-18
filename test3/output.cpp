#include <iostream>

#include "3.h"

/*void Output(const int &n, const BasicData &data,
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

struct FileName
{
    string BasicData_filename,Data_filename,Ut_file,Kt_file;
    string parameter_filename;
};
*/
void InitEnergyFile(const FileName &filename)
{
    FILE *fp = fopen(filename.Et_file.c_str(), "w");
    if (!fp) {
        printf("Cannot open Et file!\n");
        return;
    }
    fprintf(fp, "%8s %8s %16s %16s %16s %16s %16s\n",
            "time", "n", "U_all", "K_all", "F_all_x", "F_all_y", "F_all_z");
    fclose(fp);
}

void Output(const Data &data, const FileName &filename)
{
    FILE *fp1 = fopen(filename.Data_filename.c_str(), "a+");
    //n
    fprintf(fp1, "%d\n", data.n);
    //t & E
    fprintf(fp1, "   time=   %8f (fs)  Energy=  %8f (eV)\n", data.T, data.E);
    //box x*3 y*3 z*3
    fprintf(fp1, "BOX %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            data.Box.a00, data.Box.a10, data.Box.a20,  // x
            data.Box.a01, data.Box.a11, data.Box.a21,  // y
            data.Box.a02, data.Box.a12, data.Box.a22); // z
    //atoms
    for (int i = 0; i < data.n; i++)
    {
        const Atom &d = data.atoms[i];
        fprintf(fp1, " %4s  %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
                d.name.c_str(),
                d.r.a00, d.r.a10, d.r.a20,
                d.v.a00, d.v.a10, d.v.a20,
                d.dv.a00, d.dv.a10, d.dv.a20);
    }
    fclose(fp1);
    //--------------------------------------------------------------------------------------
    FILE *fp2 = fopen(filename.Et_file.c_str(), "a+");
    fprintf(fp2, "%8f %8d %16.8f %16.8f %16.8f %16.8f %16.8f\n",
             data.T,data.n,data.U_all,data.K_all,data.F_all.a00,data.F_all.a10,data.F_all.a20);
    fclose(fp2);
};
