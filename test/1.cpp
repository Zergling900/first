#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// ----------------------------
struct BasicData
{
    int n;          // number of atoms
    double T;       // unknown
    double E;       // energy
    double Box_L;   // box size
    double a0;      // 
    double theta;   // 
    double phi;     // 
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
    data.phi   *= deg2rad;
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
// R(theta, phi)
Matrix3d make_rotation(double theta, double phi)
{
    double ct = std::cos(theta);
    double st = std::sin(theta);
    double cp = std::cos(phi);
    double sp = std::sin(phi);

    Matrix3d Rz;
    Rz <<  cp, -sp, 0,
           sp,  cp, 0,
           0 ,   0 , 1;

    Matrix3d Ry;
    Ry <<  ct, 0, st,
           0 , 1, 0,
          -st, 0, ct;

    Matrix3d R = Ry * Rz;
    return R;
}


// ----------------------------
// 用一个晶格方向向量 v (比如 A_rot.col(0))
// 来估计“往该方向走一格, 在任意坐标轴上
// 会膨胀出去多少”。
// 我们用它来估算循环上限 nx,ny,nz
// ----------------------------
double span_per_step(const Vector3d &v)
{
    double sx = std::abs(v(0));
    double sy = std::abs(v(1));
    double sz = std::abs(v(2));
    return std::max({sx, sy, sz});
}


// ----------------------------
// 根据 A_rot 和 Box_L 来估计 nx,ny,nz 的上限
// 我们多加一点余量 (+2) ，然后在生成时再用 if 筛选
// ----------------------------
void estimate_cells(const Matrix3d &A_rot,
                    double Box_L,
                    int &nx_est, int &ny_est, int &nz_est)
{
    Vector3d a1 = A_rot.col(0);
    Vector3d a2 = A_rot.col(1);
    Vector3d a3 = A_rot.col(2);

    double dx = span_per_step(a1);
    double dy = span_per_step(a2);
    double dz = span_per_step(a3);

    if (dx < 1e-12) dx = 1e-12;
    if (dy < 1e-12) dy = 1e-12;
    if (dz < 1e-12) dz = 1e-12;

    nx_est = static_cast<int>( std::ceil(Box_L / dx) ) + 2;
    ny_est = static_cast<int>( std::ceil(Box_L / dy) ) + 2;
    nz_est = static_cast<int>( std::ceil(Box_L / dz) ) + 2;
}

// calculation method
//   1. A_rot = R * A
//   2. calculate nx,ny,nz
//   3. for * 3
//   4. r = A_rot * ( (ix,iy,iz) + basis[b] ) + r0
//   5. if(in Box_L) -> push_back


// ------------------------------------------------------------
// bcc
// 2
// (0,0,0) 
//(1/2,1/2,1/2)
// ------------------------------------------------------------
void bcc(FirstMolecularData &MolecularData, BasicData &data, BCCData &Data_dummy)
{
    // 随机起点 r0
    Vector3d r0(MolecularData.x0, MolecularData.y0, MolecularData.z0);

    // 晶格矩阵(立方, 未旋转)
    Matrix3d A;
    A << data.a0, 0,        0,
         0,       data.a0,  0,
         0,       0,        data.a0;

    // 旋转矩阵
    Matrix3d Rmat = make_rotation(data.theta, data.phi);

    // 旋转后的晶格矩阵
    Matrix3d A_rot = Rmat * A;

    // BCC 基底 (两个点)
    struct BCCBasis { double fx, fy, fz; };
    BCCBasis basis[2] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };

    // 粗略估计我们需要的循环上限
    int nx_est, ny_est, nz_est;
    estimate_cells(A_rot, data.Box_L, nx_est, ny_est, nz_est);

    // 清空全局 BCC，并预留空间
    BCC.clear();
    BCC.reserve(nx_est * ny_est * nz_est * 2);

    // 三重循环生成
    for (int ix = 0; ix < nx_est; ++ix)
    {
        for (int iy = 0; iy < ny_est; ++iy)
        {
            for (int iz = 0; iz < nz_est; ++iz)
            {
                Vector3d n_cell(ix, iy, iz);

                for (int b = 0; b < 2; ++b)
                {
                    Vector3d frac(n_cell(0) + basis[b].fx,
                                  n_cell(1) + basis[b].fy,
                                  n_cell(2) + basis[b].fz);

                    Vector3d r_cart = A_rot * frac + r0;

                    double x = r_cart(0);
                    double y = r_cart(1);
                    double z = r_cart(2);

                    // 裁剪到盒子里
                    if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                        x < data.Box_L &&
                        y < data.Box_L &&
                        z < data.Box_L)
                    {
                        BCCData atom;
                        atom.x1 = x;
                        atom.y1 = y;
                        atom.z1 = z;

                        atom.xv1 = 0.0;
                        atom.yv1 = 0.0;
                        atom.zv1 = 0.0;

                        atom.xdv1 = 0.0;
                        atom.ydv1 = 0.0;
                        atom.zdv1 = 0.0;

                        BCC.push_back(atom);
                    }
                }
            }
        }
    }

    // 最终原子数写回 data.n (或者你想放别处也行)
    data.n = static_cast<int>(BCC.size());

    /*  // 输出到文件的例子 (整段保留注释块)
    {
        std::ofstream fout("bcc_out.xyz");
        fout << BCC.size() << "\n";
        fout << "BCC structure\n";
        for (size_t i = 0; i < BCC.size(); ++i)
        {
            fout << "W "
                 << BCC[i].x1 << " "
                 << BCC[i].y1 << " "
                 << BCC[i].z1 << "\n";
        }
    }
    */
}


// ------------------------------------------------------------
// fcc
// 4
// (0,0,0)
// (1/2,1/2,0)
// (1/2,0,1/2)
// (0,1/2,1/2)
// ------------------------------------------------------------
void fcc(FirstMolecularData &MolecularData, BasicData &data, FCCData &Data_dummy)
{
    Vector3d r0(MolecularData.x0, MolecularData.y0, MolecularData.z0);

    Matrix3d A;
    A << data.a0, 0,        0,
         0,       data.a0,  0,
         0,       0,        data.a0;

    Matrix3d Rmat = make_rotation(data.theta, data.phi);
    Matrix3d A_rot = Rmat * A;

    struct FCCBasis { double fx, fy, fz; };
    FCCBasis basis[4] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}
    };

    int nx_est, ny_est, nz_est;
    estimate_cells(A_rot, data.Box_L, nx_est, ny_est, nz_est);

    FCC.clear();
    FCC.reserve(nx_est * ny_est * nz_est * 4);

    for (int ix = 0; ix < nx_est; ++ix)
    {
        for (int iy = 0; iy < ny_est; ++iy)
        {
            for (int iz = 0; iz < nz_est; ++iz)
            {
                Vector3d n_cell(ix, iy, iz);

                for (int b = 0; b < 4; ++b)
                {
                    Vector3d frac(n_cell(0) + basis[b].fx,
                                  n_cell(1) + basis[b].fy,
                                  n_cell(2) + basis[b].fz);

                    Vector3d r_cart = A_rot * frac + r0;

                    double x = r_cart(0);
                    double y = r_cart(1);
                    double z = r_cart(2);

                    if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                        x < data.Box_L &&
                        y < data.Box_L &&
                        z < data.Box_L)
                    {
                        FCCData atom;
                        atom.x2 = x;
                        atom.y2 = y;
                        atom.z2 = z;

                        atom.xv2 = 0.0;
                        atom.yv2 = 0.0;
                        atom.zv2 = 0.0;

                        atom.xdv2 = 0.0;
                        atom.ydv2 = 0.0;
                        atom.zdv2 = 0.0;

                        FCC.push_back(atom);
                    }
                }
            }
        }
    }

    data.n = static_cast<int>(FCC.size());

    /*  // 输出到文件的例子
    {
        std::ofstream fout("fcc_out.xyz");
        fout << FCC.size() << "\n";
        fout << "FCC structure\n";
        for (size_t i = 0; i < FCC.size(); ++i)
        {
            fout << "W "
                 << FCC[i].x2 << " "
                 << FCC[i].y2 << " "
                 << FCC[i].z2 << "\n";
        }
    }
    */
}


// ------------------------------------------------------------
// diamond
// 8
//  (0,0,0)
//  (1/2,1/2,0)
//  (1/2,0,1/2)
//  (0,1/2,1/2)
//  (1/4,1/4,1/4)
//  (3/4,3/4,1/4)
//  (3/4,1/4,3/4)
//  (1/4,3/4,3/4)
// ------------------------------------------------------------
void diamond(FirstMolecularData &MolecularData, BasicData &data, DiamondData &Data_dummy)
{
    Vector3d r0(MolecularData.x0, MolecularData.y0, MolecularData.z0);

    Matrix3d A;
    A << data.a0, 0,        0,
         0,       data.a0,  0,
         0,       0,        data.a0;

    Matrix3d Rmat = make_rotation(data.theta, data.phi);
    Matrix3d A_rot = Rmat * A;

    struct DiamondBasis { double fx, fy, fz; };
    DiamondBasis basis[8] = {
        {0.0 , 0.0 , 0.0},
        {0.5 , 0.5 , 0.0},
        {0.5 , 0.0 , 0.5},
        {0.0 , 0.5 , 0.5},
        {0.25, 0.25, 0.25},
        {0.75, 0.75, 0.25},
        {0.75, 0.25, 0.75},
        {0.25, 0.75, 0.75}
    };

    int nx_est, ny_est, nz_est;
    estimate_cells(A_rot, data.Box_L, nx_est, ny_est, nz_est);

    Diamond.clear();
    Diamond.reserve(nx_est * ny_est * nz_est * 8);

    for (int ix = 0; ix < nx_est; ++ix)
    {
        for (int iy = 0; iy < ny_est; ++iy)
        {
            for (int iz = 0; iz < nz_est; ++iz)
            {
                Vector3d n_cell(ix, iy, iz);

                for (int b = 0; b < 8; ++b)
                {
                    Vector3d frac(n_cell(0) + basis[b].fx,
                                  n_cell(1) + basis[b].fy,
                                  n_cell(2) + basis[b].fz);

                    Vector3d r_cart = A_rot * frac + r0;

                    double x = r_cart(0);
                    double y = r_cart(1);
                    double z = r_cart(2);

                    if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                        x < data.Box_L &&
                        y < data.Box_L &&
                        z < data.Box_L)
                    {
                        DiamondData atom;
                        atom.x3 = x;
                        atom.y3 = y;
                        atom.z3 = z;

                        atom.xv3 = 0.0;
                        atom.yv3 = 0.0;
                        atom.zv3 = 0.0;

                        atom.xdv3 = 0.0;
                        atom.ydv3 = 0.0;
                        atom.zdv3 = 0.0;

                        Diamond.push_back(atom);
                    }
                }
            }
        }
    }

    data.n = static_cast<int>(Diamond.size());

    /*  // 输出到文件的例子
    {
        std::ofstream fout("diamond_out.xyz");
        fout << Diamond.size() << "\n";
        fout << "Diamond structure\n";
        for (size_t i = 0; i < Diamond.size(); ++i)
        {
            fout << "C "
                 << Diamond[i].x3 << " "
                 << Diamond[i].y3 << " "
                 << Diamond[i].z3 << "\n";
        }
    }
    */
}


// ------------------------------------------------------------
// main program
// ------------------------------------------------------------
int main()
{
    BasicData data;
    read(data);

    // 给分子起始点一个随机的小平移
    FirstMolecularData MolecularData;
    random(MolecularData, data.a0);

    cout << "You have input: "
         << "E = "      << data.E
         << ", a0 = "   << data.a0
         << ", Box_L = "<< data.Box_L
         << ", theta = "<< data.theta
         << ", phi = "  << data.phi
         << std::endl;

    // 生成 BCC / FCC / Diamond 任选其一
    // 你可以根据需要调用其中一个
    {
        BCCData dummy_bcc; // 只是为了保持你的函数签名
        bcc(MolecularData, data, dummy_bcc);

        cout << "BCC atoms generated: " << BCC.size() << std::endl;
    }

    {
        FCCData dummy_fcc;
        fcc(MolecularData, data, dummy_fcc);

        cout << "FCC atoms generated: " << FCC.size() << std::endl;
    }

    {
        DiamondData dummy_dia;
        diamond(MolecularData, data, dummy_dia);

        cout << "Diamond atoms generated: " << Diamond.size() << std::endl;
    }

    cout << "Final data.n = " << data.n << std::endl;

    return 0;
}
