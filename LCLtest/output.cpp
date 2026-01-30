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
    fprintf(fp, "   time=   %f (fs)  Energy=  %f (eV)\n", data.t, data.E);
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
    fprintf(fp, "%8s %8s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n",
            "time(fs)", "n", "E", "H", "T", "s","U_all", "K_all", "F_all_x", "F_all_y", "F_all_z", "P_all_x", "P_all_y", "P_all_z");
    fclose(fp);
}

void OutputData(const Data &data, const FileName &filename, const parameter1 &pr1)
{   
    double mw1 = 1.0 / pr1.mw;
    double mb1 = 1.0 / pr1.mb;
    double s0 = data.s0;

    double E = data.E;
    double E_0 = pr1.E0;


    FILE *fp1 = fopen(filename.Data_filename.c_str(), "a+");
    //n
    fprintf(fp1, "%d\n", data.n);
    //t & E
    E = E * E_0; //to eV
    fprintf(fp1, "   time=   %8f (fs)  Energy=  %8f (eV)\n", data.t, E);
    Matrix33 Box = data.Box;

    //box x*3 y*3 z*3
    fprintf(fp1, "BOX %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            Box.a00, Box.a10, Box.a20,  // x
            Box.a01, Box.a11, Box.a21,  // y
            Box.a02, Box.a12, Box.a22); // z
    //atoms
    Matrix31 v, dv;
    for (int i = 0; i < data.n; i++)
    {
        string name = data.atoms[i].name;
        Matrix31 r = data.atoms[i].r;
        Matrix31 p = data.atoms[i].p;
        Matrix31 f = data.atoms[i].f;


         
        if (name == "W"){
            v = p * mw1 *(1.0/ s0);
            dv = f * mw1;
        }
        else if (name == "Be"){
            v = p * mb1 *(1.0/ s0);
            dv = f * mb1;
        }

        fprintf(fp1, " %4s  %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
                name.c_str(), r.a00, r.a10, r.a20, v.a00, v.a10, v.a20, dv.a00, dv.a10, dv.a20);
    }
    fclose(fp1);
    //--------------------------------------------------------------------------------------
}

void OutputEnergy(const Data &data, const FileName &filename, const parameter1 &pr1)
{
    FILE *fp2 = fopen(filename.Et_file.c_str(), "a+");
    int N = data.n;
    double E = data.E;
    double U = data.U_all;
    double K = data.K_all;
    double H = data.H;
    double T = data.T;

    double F_x = data.F_all.a00;
    double F_y = data.F_all.a10;
    double F_z = data.F_all.a20;

    double P_x = data.P_all.a00;
    double P_y = data.P_all.a10;
    double P_z = data.P_all.a20;

    double s0 = data.s0;
    double t = data.t; // fs
    double E_0 = pr1.E0;


    E = E *E_0; //to eV
    U = U *E_0;
    K = K *E_0;
    H = H *E_0;

    fprintf(fp2, "%8f %8d %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
             t,N, E, H, T, s0, U,K,F_x,F_y,F_z,P_x,P_y,P_z);
    fclose(fp2);
};
