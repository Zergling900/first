#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"
struct Cell_List
{
    double Wx,Wy,Wz;
    int Mx,My,Mz,cell_num;
    std::vector<int> Cell;//N
    //std::vector<int> num_in_cell;  //M number of particles in each cell
    std::vector<int> cell_offset;  // The starting point of cell[i] in [atom_indices]
    std::vector<int> atom_indices; // Install the particle IDs in segments into cells(每个粒子在所在的格子中的位置（第几个）)
    //std::vector<int> atom_order;   // Local order of particle[i] within its cell
};
void lcl0(Data&data,Cell_List&cl,const parameter1 &pr1, const parameter2 &pr2)
{
    cl.Mx = floor(data.Box.a00 / pr2.R);
    cl.My = floor(data.Box.a11 / pr2.R);
    cl.Mz = floor(data.Box.a22 / pr2.R);

    cl.Wx = (data.Box.a00/cl.Mx);
    cl.Wy = (data.Box.a11/cl.My);
    cl.Wz = (data.Box.a22/cl.Mz);

    cl.cell_num = cl.Mx * cl.My * cl.Mz;

    cl.Cell.resize(data.n);
    cl.cell_offset.resize(cl.cell_num + 1);
    cl.atom_indices.resize(data.n);
}

void lcl1(Data&data,Cell_List&cl,const parameter1 &pr1, const parameter2 &pr2)
{
    double Wx = cl.Wx;
    double Wy = cl.Wy;
    double Wz = cl.Wz;
    std::vector<int> num_in_cell(cl.cell_num, 0);

    for(int i=0;i<data.n;i++)
    {
        double mx = floor(data.atoms[i].r.a00 / Wx);
        double my = floor(data.atoms[i].r.a10 / Wy);
        double mz = floor(data.atoms[i].r.a20 / Wz);
        cl.Cell[i] = mx + my * cl.Mx + mz * cl.Mx * cl.My;
        num_in_cell[cl.Cell[i]] ++; 
    }
    
    cl.cell_offset[0] = 0;
    for(int i=0;i<cl.cell_num;i++)
    {
        cl.cell_offset[i+1] = cl.cell_offset[i] + num_in_cell[i];
    }

    std::vector<int> offset1(cl.cell_num);
    for(int i=0;i<cl.cell_num;i++)
    {
        offset1[i] = cl.cell_offset[i];
    }

    for(int i=0;i<data.n;i++)
    {
        int c = cl.Cell[i];//读取第i个粒子所在的格子编号
        int pos = offset1[c]++;//该格子中粒子的个数加1，并返回该粒子在该格子中的位置
        cl.atom_indices[pos] = i;//将粒子i的编号存入该格子中
    }



}