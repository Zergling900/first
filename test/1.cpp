#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// ----------------------------
struct BasicData
{
    int n;        // number of atoms
    double T;     // unknown
    double E;     // energy
    double Box_L; // box size
    double a0;    //
    double theta; //
    double phi;   //
};

struct FirstMolecularData
{
    double x0, y0, z0;
};

// BCCData
struct BCCData
{
    double x1, y1, z1;
    double xv1, yv1, zv1;
    double xdv1, ydv1, zdv1;
};
std::vector<BCCData> BCC;
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

// FCCData
struct FCCData
{
    double x2, y2, z2;
    double xv2, yv2, zv2;
    double xdv2, ydv2, zdv2;
};
std::vector<FCCData> FCC;

// DiamondData
struct DiamondData
{
    double x3, y3, z3;
    double xv3, yv3, zv3;
    double xdv3, ydv3, zdv3;
};
std::vector<DiamondData> Diamond;

// ----------------------------

void read(BasicData &data)
{
    cout << "please input energy 'E', cubic side length 'a0', box size 'Box_L', 'theta', 'phi': " << std::endl;
    cin >> data.E >> data.a0 >> data.Box_L >> data.theta >> data.phi;
    // ----------------
    data.n = 0;
    data.T = 0.0;

    // angle to radian
    const double deg2rad = M_PI / 180.0;
    data.theta *= deg2rad;
    data.phi *= deg2rad;
}

// ----------------------------
// random (first molecular)
std::mt19937 global_eng(std::random_device{}());
void random(FirstMolecularData &MolecularData, double a0)
{
    std::uniform_real_distribution<> distr(0.0, 0.1 * a0);
    MolecularData.x0 = distr(global_eng);
    MolecularData.y0 = distr(global_eng);
    MolecularData.z0 = distr(global_eng);
}

// ----------------------------

// Rotation matrix (theta, phi)
Matrix3d Rotation(double theta, double phi)
{
    double ct = std::cos(theta);
    double st = std::sin(theta);
    double cp = std::cos(phi);
    double sp = std::sin(phi);

    Matrix3d Rz;
    Rz << cp, -sp, 0,
        sp, cp, 0,
        0, 0, 1;

    Matrix3d Ry;
    Ry << ct, 0, st,
        0, 1, 0,
        -st, 0, ct;

    Matrix3d R = Ry * Rz;
    return R;
}

// basic matrix A
Matrix3d A0(const BasicData &data)
{
    Matrix3d A;
    A << data.a0, 0, 0,
        0, data.a0, 0,
        0, 0, data.a0;
    return A;
}

// ----------------------------

double Step(const Vector3d &a)
{
    double ax = std::abs(a(0));
    double ay = std::abs(a(1));
    double az = std::abs(a(2));
    return std::max({ax, ay, az});
}

// ----------------------------
/*
a0 = a0 0 0
     0  a0 0
     0  0  a0

A_rot = R * A = Ry * Rz *A
      =  ct, 0, st,        cp, -sp, 0         a0, 0, 0
         0 , 1, 0,    X    sp,  cp, 0    X    0, a0, 0
        -st, 0, ct;        0 ,  0,  1.        0, 0, a0
*/
// ----------------------------
void Box(const Matrix3d &A_rot, const BasicData &data, int &nx_est, int &ny_est, int &nz_est)
{
    // col(0) = (a0x, a0y, a0z)
    Vector3d a1 = A_rot.col(0);
    Vector3d a2 = A_rot.col(1);
    Vector3d a3 = A_rot.col(2);

    double dx = Step(a1);
    double dy = Step(a2);
    double dz = Step(a3);

    if (dx < 1e-6)
        dx = 1e-6;
    if (dy < 1e-6)
        dy = 1e-6;
    if (dz < 1e-6)
        dz = 1e-6;

    nx_est = static_cast<int>(std::ceil(data.Box_L / dx)) + 1;
    ny_est = static_cast<int>(std::ceil(data.Box_L / dy)) + 1;
    nz_est = static_cast<int>(std::ceil(data.Box_L / dz)) + 1;
}
// calculation method
//   1. A_rot = R * A
//   2. calculate nx,ny,nz
//   3. for * 3
//   4. r = A_rot * ( (ix,iy,iz) + basis[b] ) + r0
//   5. if(in Box_L) -> push_back

// ------------------------------------------------------------
void bcc(FirstMolecularData &MolecularData, BasicData &data, BCCData &Data)
{
    // random first atom position r0
    Vector3d r0(MolecularData.x0, MolecularData.y0, MolecularData.z0);

    // rotation matrix
    Matrix3d R = Rotation(data.theta, data.phi);
    Matrix3d A = A0(data);
    Matrix3d A_rot = R * A;
    // BCC basis
    // 2
    // (0,0,0)
    //(1/2,1/2,1/2)
    struct BCCBasis
    {
        double fx, fy, fz;
    };
    BCCBasis basis[2] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}};

    // Box calculation
    int nx_est, ny_est, nz_est;
    Box(A_rot, data, nx_est, ny_est, nz_est);

    // one bcc cell has 2 atoms
    BCC.clear();
    BCC.reserve(nx_est * ny_est * nz_est * 2);

    // --------------------------------------------
    for (int ix = 0; ix <= nx_est; ++ix)
    {
        for (int iy = 0; iy <= ny_est; ++iy)
        {
            for (int iz = 0; iz <= nz_est; ++iz)
            {
                Vector3d i(ix, iy, iz);

                for (int b = 0; b < 2; ++b) // bcc 2 basis atoms
                {
                    Vector3d frac(i(0) + basis[b].fx, // fractional coordinates
                                  i(1) + basis[b].fy,
                                  i(2) + basis[b].fz);

                    Vector3d r_cart = A_rot * frac + r0;

                    double x = r_cart(0);
                    double y = r_cart(1);
                    double z = r_cart(2);

                    // in box
                    if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                        x < data.Box_L &&
                        y < data.Box_L &&
                        z < data.Box_L)
                    {
                        BCCData bcc;
                        bcc.x1 = x;
                        bcc.y1 = y;
                        bcc.z1 = z;

                        // bcc.xv1 = bcc.xdv1;
                        // bcc.yv1 = bcc.ydv1;
                        // bcc.zv1 = bcc.zdv1;
                        bcc.xv1 = 0.0;
                        bcc.yv1 = 0.0;
                        bcc.zv1 = 0.0;

                        bcc.xdv1 = 0.0;
                        bcc.ydv1 = 0.0;
                        bcc.zdv1 = 0.0;

                        BCC.push_back(bcc);
                    }
                }
            }
        }
    }

    data.n = static_cast<int>(BCC.size());
}

void output_bcc(const BasicData &data, const std::vector<BCCData> &BCC,
                const std::string &filename = "bcc.md3")
{
    std::ofstream fout(filename);
    fout << data.n << "\n";
    fout << "   time=   " << data.T
         << " (fs)  Energy=  " << data.E
         << " (eV)" << "\n";

    fout << "BOX"
         << std::fixed << std::setprecision(8)
         << std::setw(18) << data.Box_L // ax
         << std::setw(16) << 0.0        // ay
         << std::setw(16) << 0.0        // az
         << std::setw(16) << 0.0        // bx
         << std::setw(16) << data.Box_L // by
         << std::setw(16) << 0.0        // bz
         << std::setw(16) << 0.0        // cx
         << std::setw(16) << 0.0        // cy
         << std::setw(16) << data.Box_L // cz
         << "\n";

    for (size_t i = 0; i < data.n; ++i)
    {
        fout << std::setw(4) << "W"
             << std::fixed << std::setprecision(8)
             << std::setw(18) << BCC[i].x1
             << std::setw(16) << BCC[i].y1
             << std::setw(16) << BCC[i].z1
             << std::setw(16) << BCC[i].xv1
             << std::setw(16) << BCC[i].yv1
             << std::setw(16) << BCC[i].zv1
             << std::setw(16) << BCC[i].xdv1
             << std::setw(16) << BCC[i].ydv1
             << std::setw(16) << BCC[i].zdv1
             << "\n";
    }

    fout.close();
}
// ------------------------------------------------------------
// fcc
// 4
// (0,0,0)
// (1/2,1/2,0)
// (1/2,0,1/2)
// (0,1/2,1/2)
// ------------------------------------------------------------
void fcc(FirstMolecularData &MolecularData, BasicData &data, FCCData &Data)
{
    // random first atom position r0
    Vector3d r0(MolecularData.x0, MolecularData.y0, MolecularData.z0);

    // rotation matrix
    Matrix3d R = Rotation(data.theta, data.phi);
    Matrix3d A = A0(data);
    Matrix3d A_rot = R * A;
    // FCC basis
    // 4
    // (0,0,0)
    // (1/2,1/2,0)
    // (1/2,0,1/2)
    // (0,1/2,1/2)
    struct FCCBasis
    {
        double fx, fy, fz;
    };
    FCCBasis basis[4] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}};

    // Box calculation
    int nx_est, ny_est, nz_est;
    Box(A_rot, data, nx_est, ny_est, nz_est);

    // one fcc cell has 4 atoms
    FCC.clear();
    FCC.reserve(nx_est * ny_est * nz_est * 4);

    // --------------------------------------------
    for (int ix = 0; ix <= nx_est; ++ix)
    {
        for (int iy = 0; iy <= ny_est; ++iy)
        {
            for (int iz = 0; iz <= nz_est; ++iz)
            {
                Vector3d i(ix, iy, iz);

                for (int b = 0; b < 4; ++b) // fcc 4 basis atoms
                {
                    Vector3d frac(i(0) + basis[b].fx, // fractional coordinates
                                  i(1) + basis[b].fy,
                                  i(2) + basis[b].fz);

                    Vector3d r_cart = A_rot * frac + r0;

                    double x = r_cart(0);
                    double y = r_cart(1);
                    double z = r_cart(2);

                    // in box
                    if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                        x < data.Box_L &&
                        y < data.Box_L &&
                        z < data.Box_L)
                    {
                        FCCData fcc;
                        fcc.x2 = x;
                        fcc.y2 = y;
                        fcc.z2 = z;

                        // fcc.xv2 = fcc.xdv2;
                        // fcc.yv2 = fcc.ydv2;
                        // fcc.zv2 = fcc.zdv2;
                        fcc.xv2 = 0.0;
                        fcc.yv2 = 0.0;
                        fcc.zv2 = 0.0;

                        fcc.xdv2 = 0.0;
                        fcc.ydv2 = 0.0;
                        fcc.zdv2 = 0.0;

                        FCC.push_back(fcc);
                    }
                }
            }
        }
    }

    data.n = static_cast<int>(FCC.size());
}
void output_fcc(const BasicData &data, const std::vector<FCCData> &FCC,
                const std::string &filename = "fcc.md3")
{
    std::ofstream fout(filename);
    fout << data.n << "\n";
    fout << "   time=   " << data.T
         << " (fs)  Energy=  " << data.E
         << " (eV)" << "\n";

    fout << "BOX"
         << std::fixed << std::setprecision(8)
         << std::setw(18) << data.Box_L // ax
         << std::setw(16) << 0.0        // ay
         << std::setw(16) << 0.0        // az
         << std::setw(16) << 0.0        // bx
         << std::setw(16) << data.Box_L // by
         << std::setw(16) << 0.0        // bz
         << std::setw(16) << 0.0        // cx
         << std::setw(16) << 0.0        // cy
         << std::setw(16) << data.Box_L // cz
         << "\n";

    for (size_t i = 0; i < data.n; ++i)
    {
        fout << std::setw(4) << "Cu"
             << std::fixed << std::setprecision(8)
             << std::setw(18) << FCC[i].x2
             << std::setw(16) << FCC[i].y2
             << std::setw(16) << FCC[i].z2
             << std::setw(16) << FCC[i].xv2
             << std::setw(16) << FCC[i].yv2
             << std::setw(16) << FCC[i].zv2
             << std::setw(16) << FCC[i].xdv2
             << std::setw(16) << FCC[i].ydv2
             << std::setw(16) << FCC[i].zdv2
             << "\n";
    }

    fout.close();
}
// ------------------------------------------------------------
// diamond
// 8
// (0,0,0)
// (1/2,1/2,0)
// (1/2,0,1/2)
// (0,1/2,1/2)
// (1/4,1/4,1/4)
// (3/4,3/4,1/4)
// (3/4,1/4,3/4)
// (1/4,3/4,3/4)
// ------------------------------------------------------------

void diamond(FirstMolecularData &MolecularData, BasicData &data, DiamondData &Data)
{
    // random first atom position r0
    Vector3d r0(MolecularData.x0, MolecularData.y0, MolecularData.z0);

    // rotation matrix
    Matrix3d R = Rotation(data.theta, data.phi);
    Matrix3d A = A0(data);
    Matrix3d A_rot = R * A;
    // diamond
    // 8
    // (0,0,0)
    // (1/2,1/2,0)
    // (1/2,0,1/2)
    // (0,1/2,1/2)
    // (1/4,1/4,1/4)
    // (3/4,3/4,1/4)
    // (3/4,1/4,3/4)
    // (1/4,3/4,3/4)
    struct DiamondBasis
    {
        double fx, fy, fz;
    };
    DiamondBasis basis[8] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5},
        {0.25, 0.25, 0.25},
        {0.75, 0.75, 0.25},
        {0.75, 0.25, 0.75},
        {0.25, 0.75, 0.75}};

    // Box calculation
    int nx_est, ny_est, nz_est;
    Box(A_rot, data, nx_est, ny_est, nz_est);

    // one diamond cell has 8 atoms
    Diamond.clear();
    Diamond.reserve(nx_est * ny_est * nz_est * 8);

    // --------------------------------------------
    for (int ix = 0; ix <= nx_est; ++ix)
    {
        for (int iy = 0; iy <= ny_est; ++iy)
        {
            for (int iz = 0; iz <= nz_est; ++iz)
            {
                Vector3d i(ix, iy, iz);

                for (int b = 0; b < 8; ++b) // diamond 8 basis atoms
                {
                    Vector3d frac(i(0) + basis[b].fx, // fractional coordinates
                                  i(1) + basis[b].fy,
                                  i(2) + basis[b].fz);

                    Vector3d r_cart = A_rot * frac + r0;

                    double x = r_cart(0);
                    double y = r_cart(1);
                    double z = r_cart(2);

                    // in box
                    if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                        x < data.Box_L &&
                        y < data.Box_L &&
                        z < data.Box_L)
                    {
                        DiamondData dia;
                        dia.x3 = x;
                        dia.y3 = y;
                        dia.z3 = z;

                        // dia.xv3 = dia.xdv3;
                        // dia.yv3 = dia.ydv3;
                        // dia.zv3 = dia.zdv3;
                        dia.xv3 = 0.0;
                        dia.yv3 = 0.0;
                        dia.zv3 = 0.0;

                        dia.xdv3 = 0.0;
                        dia.ydv3 = 0.0;
                        dia.zdv3 = 0.0;

                        Diamond.push_back(dia);
                    }
                }
            }
        }
    }

    data.n = static_cast<int>(Diamond.size());
}
void output_diamond(const BasicData &data, const std::vector<DiamondData> &Diamond,
                const std::string &filename = "diamond.md3")
{
    std::ofstream fout(filename);
    fout << data.n << "\n";
    fout << "   time=   " << data.T
         << " (fs)  Energy=  " << data.E
         << " (eV)" << "\n";

    fout << "BOX"
         << std::fixed << std::setprecision(8)
         << std::setw(18) << data.Box_L // ax
         << std::setw(16) << 0.0        // ay
         << std::setw(16) << 0.0        // az
         << std::setw(16) << 0.0        // bx
         << std::setw(16) << data.Box_L // by
         << std::setw(16) << 0.0        // bz
         << std::setw(16) << 0.0        // cx
         << std::setw(16) << 0.0        // cy
         << std::setw(16) << data.Box_L // cz
         << "\n";

    for (size_t i = 0; i < data.n; ++i)
    {
        fout << std::setw(4) << "C"
             << std::fixed << std::setprecision(8)
             << std::setw(18) << Diamond[i].x3
             << std::setw(16) << Diamond[i].y3
             << std::setw(16) << Diamond[i].z3
             << std::setw(16) << Diamond[i].xv3
             << std::setw(16) << Diamond[i].yv3
             << std::setw(16) << Diamond[i].zv3
             << std::setw(16) << Diamond[i].xdv3
             << std::setw(16) << Diamond[i].ydv3
             << std::setw(16) << Diamond[i].zdv3
             << "\n";
    }

    fout.close();
}

// ------------------------------------------------------------
// main program
// ------------------------------------------------------------
int main()
{
    BasicData data;
    FirstMolecularData firstMol;

    // 1. read
    read(data);

    // 2. first atom random
    random(firstMol, data.a0);

    // 3. bcc or fcc or diamond
    cout << "choose structure type: 1 = BCC (W), 2 = FCC (Cu), 3 = Diamond (C)" << endl;
    int choice;
    cin >> choice;

    // 4. output
    if (choice == 1)
    {
        BCCData dummy;
        bcc(firstMol, data, dummy);
        output_bcc(data, BCC, "bcc.md3");
        cout << "[BCC] generated " << data.n << "written to bcc.md3\n";
    }
    else if (choice == 2)
    {
        FCCData dummy;
        fcc(firstMol, data, dummy);
        output_fcc(data, FCC, "fcc.md3");
        cout << "[FCC] generated " << data.n << "written to fcc.md3\n";
    }
    else if (choice == 3)
    {
        DiamondData dummy;
        diamond(firstMol, data, dummy);
        output_diamond(data, Diamond, "diamond.md3");
        cout << "[Diamond] generated " << data.n << "written to diamond.md3\n";
    }
    else
    {
        cerr << "invalid choice.\n";
        return 1;
    }

    return 0;
}
