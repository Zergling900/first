#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
//
struct BasicData
{
    int n;
    double T, E, Box_L, a0, theta, phi;
};
//
struct FirstMolecularData
{
    double x0, y0, z0;
};
// BCCData
struct BCCData
{
    double x1, y1, z1, xv1, yv1, zv1, xdv1, ydv1, zdv1;
};
std::vector<BCCData> BCC;
// FCCData
struct FCCData
{
    double x2, y2, z2, xv2, yv2, zv2, xdv2, ydv2, zdv2;
};
std::vector<FCCData> FCC;
// DiamondData
struct DiamondData
{
    double x3, y3, z3, xv3, yv3, zv3, xdv3, ydv3, zdv3;
};
std::vector<DiamondData> Diamond;

// read input data(now i don't konw energy 'E')
void read(BasicData &data)
{
    cout << "please input energy 'E', cubic side length 'a0', box size 'Box_L', 'theta', 'phi': "
         << std::endl;
    cin >> data.E >> data.a0 >> data.Box_L >> data.theta >> data.phi;
}
// radom first molecular
std::mt19937 global_eng(std::random_device{}());
void random(FirstMolecularData &MolecularData, double a0)
{
    std::uniform_real_distribution<> distr(0.0, 0.1 * a0);
    MolecularData.x0 = distr(global_eng);
    MolecularData.y0 = distr(global_eng);
    MolecularData.z0 = distr(global_eng);
}

/*
bcc
(0,0,0),(0,0,1),(0,1,0),(1,0,0),(1,1,0),(1,0,1),(0,1,1),(1,1,1),(1/2,1/2,1/2)
*/
void bcc(FirstMolecularData &MolecularData, BasicData &data, BCCData &Data)
{
    int nx1 = ((data.Box_L - MolecularData.x0) / data.a0);
    int ny1 = ((data.Box_L - MolecularData.y0) / data.a0);
    int nz1 = ((data.Box_L - MolecularData.z0) / data.a0);

    int n1 = nx1 * ny1 * nz1;

    vector <BCCData> BCC = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };

    for (int ix = 0; ix < nx1-1; ++ix)
    {
        for (int iy = 0; iy < ny1-1; ++iy)
        {
            for (int iz = 0; iz < nz1-1; ++iz)
            {
                {
                    BCCData bcc;
                    bcc.x1 = MolecularData.x0 + data.a0 * (ix + 0.0);
                    bcc.y1 = MolecularData.y0 + data.a0 * (iy + 0.0);
                    bcc.z1 = MolecularData.z0 + data.a0 * (iz + 0.0);

                    bcc.xv1 = 0.0;
                    bcc.yv1 = 0.0;
                    bcc.zv1 = 0.0;

                    bcc.xdv1 = 0.0;
                    bcc.ydv1 = 0.0;
                    bcc.zdv1 = 0.0;

                    BCC.push_back(bcc);
                }
                {
                    BCCData bcc;
                    bcc.x1 = MolecularData.x0 + data.a0 * (ix + 0.5);
                    bcc.y1 = MolecularData.y0 + data.a0 * (iy + 0.5);
                    bcc.z1 = MolecularData.z0 + data.a0 * (iz + 0.5);

                    bcc.xv1 = 0.0;
                    bcc.yv1 = 0.0;
                    bcc.zv1 = 0.0;

                    bcc.xdv1 = 0.0;
                    bcc.ydv1 = 0.0;
                    bcc.zdv1 = 0.0;

                    BCC.push_back(bcc);

                    if (bcc.x1 < data.Box_L &&
                        bcc.y1 < data.Box_L &&
                        bcc.z1 < data.Box_L)
                    {
                        BCC.push_back(bcc);
                    }
                }
            }
        }
    }
    n1 = BCC.size();
}

// fcc
void fcc(FirstMolecularData &MolecularData, BasicData &data, FCCData &Data)
{
    int nx2 = ((data.Box_L - MolecularData.x0) / data.a0);
    int ny2 = ((data.Box_L - MolecularData.y0) / data.a0);
    int nz2 = ((data.Box_L - MolecularData.z0) / data.a0);

    int n2 = nx2 * ny2 * nz2; // 先算一个理论值，后面我们会用 FCC.size() 覆盖掉它

    for (int ix = 0; ix < nx2; ++ix)
    {
        for (int iy = 0; iy < ny2; ++iy)
        {
            for (int iz = 0; iz < nz2; ++iz)
            {
                // -------- (0.0, 0.0, 0.0) --------
                {
                    FCCData fcc;
                    fcc.x2 = MolecularData.x0 + data.a0 * (ix + 0.0);
                    fcc.y2 = MolecularData.y0 + data.a0 * (iy + 0.0);
                    fcc.z2 = MolecularData.z0 + data.a0 * (iz + 0.0);

                    fcc.xv2 = 0.0;
                    fcc.yv2 = 0.0;
                    fcc.zv2 = 0.0;

                    fcc.xdv2 = 0.0;
                    fcc.ydv2 = 0.0;
                    fcc.zdv2 = 0.0;

                    // 边界检查，防止超出 Box_L
                    if (fcc.x2 < data.Box_L &&
                        fcc.y2 < data.Box_L &&
                        fcc.z2 < data.Box_L)
                    {
                        FCC.push_back(fcc);
                    }
                }

                // -------- (0.5, 0.5, 0.0) --------
                {
                    FCCData fcc;
                    fcc.x2 = MolecularData.x0 + data.a0 * (ix + 0.5);
                    fcc.y2 = MolecularData.y0 + data.a0 * (iy + 0.5);
                    fcc.z2 = MolecularData.z0 + data.a0 * (iz + 0.0);

                    fcc.xv2 = 0.0;
                    fcc.yv2 = 0.0;
                    fcc.zv2 = 0.0;

                    fcc.xdv2 = 0.0;
                    fcc.ydv2 = 0.0;
                    fcc.zdv2 = 0.0;

                    if (fcc.x2 < data.Box_L &&
                        fcc.y2 < data.Box_L &&
                        fcc.z2 < data.Box_L)
                    {
                        FCC.push_back(fcc);
                    }
                }
                // -------- (0.5, 0.0, 0.5) --------
                {
                    FCCData fcc;
                    fcc.x2 = MolecularData.x0 + data.a0 * (ix + 0.5);
                    fcc.y2 = MolecularData.y0 + data.a0 * (iy + 0.0);
                    fcc.z2 = MolecularData.z0 + data.a0 * (iz + 0.5);

                    fcc.xv2 = 0.0;
                    fcc.yv2 = 0.0;
                    fcc.zv2 = 0.0;

                    fcc.xdv2 = 0.0;
                    fcc.ydv2 = 0.0;
                    fcc.zdv2 = 0.0;

                    if (fcc.x2 < data.Box_L &&
                        fcc.y2 < data.Box_L &&
                        fcc.z2 < data.Box_L)
                    {
                        FCC.push_back(fcc);
                    }
                }

                // -------- (0.0, 0.5, 0.5) --------
                {
                    FCCData fcc;
                    fcc.x2 = MolecularData.x0 + data.a0 * (ix + 0.0);
                    fcc.y2 = MolecularData.y0 + data.a0 * (iy + 0.5);
                    fcc.z2 = MolecularData.z0 + data.a0 * (iz + 0.5);

                    fcc.xv2 = 0.0;
                    fcc.yv2 = 0.0;
                    fcc.zv2 = 0.0;

                    fcc.xdv2 = 0.0;
                    fcc.ydv2 = 0.0;
                    fcc.zdv2 = 0.0;

                    if (fcc.x2 < data.Box_L &&
                        fcc.y2 < data.Box_L &&
                        fcc.z2 < data.Box_L)
                    {
                        FCC.push_back(fcc);
                    }
                }
            }
        }
    }
    n2 = FCC.size();
}

// diamond
void diamond(FirstMolecularData &MolecularData, BasicData &data, DiamondData &Data)
{
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
