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
std::vector<BCCData> BCC;
//FCCData
struct FCCData
{
    double x2, y2, z2, xv2, yv2, zv2, xdv2, ydv2, zdv2;
};
std::vector<FCCData> FCC
//DiamondData
struct DiamondData
{
    double x3, y3, z3, xv3, yv3, zv3, xdv3, ydv3, zdv3;
};
std::vector<FCCData> D
// read input data(now i don't konw energy 'E')
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
    for (int ix = 0;ix < nx; ++ix)
    for(int iy = 0;iy < nx; ++iy)
    for(int iz = 0;iz < nx; ++iz){
        {
            BCC a;
            a.x = molecularData.x0 + data.a0 * (ix + 0.0);
        a.y = molecularData.y0 + data.a0 * (iy + 0.0);
        a.z = molecularData.z0 + data.a0 * (iz + 0.0);

        a.vx = vel_dist(global_eng);
        a.vy = vel_dist(global_eng);
        a.vz = vel_dist(global_eng);

        a.ax = 0.0;
        a.ay = 0.0;
        a.az = 0.0;
        }
        
        
        
    }
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
