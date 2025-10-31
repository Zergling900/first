#include <iostream>
#include <random>

using namespace std;
//
struct BasicData {
    int n;
    double T, E, Box_L, a0;
};
//
struct FirstMolecularData
{
    double x0, y0, z0;
};
//BCCData
struct BCCData
{
    double x1, y1, z1, xv1, yv1, zv1, xdv1, ydv1, zdv1;
};
//FCCData
struct FCCData
{
    double x2, y2, z2, xv2, yv2, zv2, xdv2, ydv2, zdv2;
};
//DiamondData
struct DiamondData
{
    double x3, y3, z3, xv3, yv3, zv3, xdv3, ydv3, zdv3;
};
// read input data（now i don't konw energy 'E'）
void read(BasicData &data)
{
    cout << "please input energy 'E', cubic side length 'a0', box size 'Box_L': " 
    << std::endl;
    cin >> data.E >> data.a0 >> data.Box_L;
}
//radom first molecular 
std::mt19937 global_eng(std::random_device{}());
void random(FirstMolecularData &molecularData, double a0)
{
    std::uniform_real_distribution<> distr(0.0, a0);
    molecularData.x0 = distr(global_eng);
    molecularData.y0 = distr(global_eng);
    molecularData.z0 = distr(global_eng);
}

/*
bcc
(0,0,0),(0,0,1),(0,1,0),(1,0,0),(1,1,0),(1,0,1),(0,1,1),(1,1,1),(1/2,1/2,1/2)
*/ 
void bcc(FirstMolecularData &molecularData, BasicData &data,BCCData &Data)
{
    int n = (data.Box_L/molecularData.x0)*(data.Box_L/molecularData.y0)*(data.Box_L/molecularData.z0);

}

// fcc
void fcc(FirstMolecularData &molecularData, BasicData &data,FCCData &Data)
{
    int n = (data.Box_L/molecularData.x0)*(data.Box_L/molecularData.y0)*(data.Box_L/molecularData.z0);

}
{

}
// diamond
void diamond(FirstMolecularData &molecularData, BasicData &data, DiamondData &Data)
{
    int n = (data.Box_L/molecularData.x0)*(data.Box_L/molecularData.y0)*(data.Box_L/molecularData.z0);

}

// main program
int main()
{
    BasicData data;
    read(data);

    cout << "You have input: n = " << data.n << ", T = " << data.T
         << ", E = " << data.E << ", Box_L = " << data.Box_L << std::endl;

    // Example usage of random function
    double x, y, z;
    random(x, y, z, data.Box_L);
    cout << "Random coordinates: (" << x << ", " << y << ", " << z << ")" << std::endl;

    return 0;
}
