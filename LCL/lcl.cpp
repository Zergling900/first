#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"

void lcl0(Data&data,Cell_List&cell_list,const parameter1 &pr1, const parameter2 &pr2)
{
    int N = data.n;
    const double Lx = data.Box.a00;
    const double Ly = data.Box.a11;
    const double Lz = data.Box.a22;
    double r_cut = pr2.R;
    cell_list.Mx = floor(Lx/r_cut);
    cell_list.My = floor(Ly/r_cut);
    cell_list.Mz = floor(Lz/r_cut);
    cell_list.Wx = Lx/cell_list.Mx;
    cell_list.Wy = Ly/cell_list.My;
    cell_list.Wz = Lz/cell_list.Mz;
    cell_list.cell_num = cell_list.Mx*cell_list.My*cell_list.Mz;
    int M = cell_list.cell_num;

    double Wx = cell_list.Wx;
    double Wy = cell_list.Wy;
    double Wz = cell_list.Wz;
    int Mx = cell_list.Mx;
    int My = cell_list.My;
    int Mz = cell_list.Mz;
    int cell_num = cell_list.cell_num;

    cell_list.cell.resize(N);
    cell_list.atom_order.resize(N);
    cell_list.atom_indices.resize(N);

    cell_list.num_in_cell.resize(M, -1);
    cell_list.cell_offset.resize(M + 1);

    for (int i = 0; i < cell_num; i++)
    {
        cell_list.num_in_cell[i] = -1;
    }

    for (int i = 0; i < N; i++)
    {
        int cx  = floor(data.atoms[i].r.a00/Wx);
        int cy  = floor(data.atoms[i].r.a10/Wy);
        int cz  = floor(data.atoms[i].r.a20/Wz);
        int c_id = cx + cy*Mx + cz*Mx*My;
        cell_list.cell[i] = c_id;
        cell_list.num_in_cell[c_id]++;
        cell_list.atom_order[i] = cell_list.num_in_cell[c_id];
        /*第i个粒子对应的cell_id（cell[i])
        cell_id里包含的粒子数+1
        第i个粒子在cell_id里的序号*/  
    }
    
    
};  