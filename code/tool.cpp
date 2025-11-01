#include <Eigen/Dense>
#include <vector>
#include <cmath>    // for sin, cos
#include <iostream>

using namespace Eigen;
using namespace std;

struct AtomData {
    double x,y,z;
    double vx,vy,vz;
    double ax,ay,az;
};

// 生成旋转矩阵 R = Ry(theta) * Rz(phi)
Matrix3d make_rotation(double theta, double phi) {
    double ct = cos(theta);
    double st = sin(theta);
    double cp = cos(phi);
    double sp = sin(phi);

    // Rz(phi)
    Matrix3d Rz;
    Rz <<  cp, -sp, 0,
           sp,  cp, 0,
           0 ,   0 , 1;

    // Ry(theta)
    Matrix3d Ry;
    Ry <<  ct, 0, st,
           0 , 1, 0,
          -st, 0, ct;

    // total
    return Ry * Rz;
}

// 主例子：生成一个旋转后的钻石晶格
int main() {
    // ---------------- 参数区 ----------------
    double a0    = 3.57;    // 晶格常数
    double Box_L = 20.0;    // 盒子大小
    int nx = 4, ny = 4, nz = 4; // 晶胞复制数

    // 随便指定个整体平移（也可以是0,0,0）
    Vector3d r0(0.0, 0.0, 0.0);

    // 旋转角，单位：弧度
    double theta = 0.3; // 仰角, 比如 ~17度
    double phi   = 1.0; // 方位角, 比如 ~57度

    // 晶格矩阵 A（基本是一个立方：a0 * I）
    Matrix3d A;
    A << a0, 0,   0,
         0,  a0, 0,
         0,  0,  a0;

    // 旋转矩阵 R
    Matrix3d R = make_rotation(theta, phi);

    // 旋转后的晶格矩阵
    Matrix3d A_rot = R * A;

    // 钻石晶格基底（8个原子）
    vector<Vector3d> basis = {
        {0.0 , 0.0 , 0.0},
        {0.5 , 0.5 , 0.0},
        {0.5 , 0.0 , 0.5},
        {0.0 , 0.5 , 0.5},
        {0.25, 0.25, 0.25},
        {0.75, 0.75, 0.25},
        {0.75, 0.25, 0.75},
        {0.25, 0.75, 0.75}
    };

    vector<AtomData> atoms;
    atoms.reserve(nx * ny * nz * basis.size());

    // 生成 + 边界过滤，一步完成（高效版本）
    for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
    for (int iz = 0; iz < nz; ++iz)
    {
        Vector3d n_cell(ix, iy, iz);

        for (const auto &rb : basis) {
            // r = A_rot * (n + rb) + r0
            Vector3d frac = n_cell + rb;   // 晶胞内的分数+整数位置
            Vector3d r_cart = A_rot * frac + r0;

            double x = r_cart(0);
            double y = r_cart(1);
            double z = r_cart(2);

            // 剪在盒子里（跟你原来 if(...) 一样）
            if (x >= 0.0 && y >= 0.0 && z >= 0.0 &&
                x < Box_L && y < Box_L && z < Box_L)
            {
                AtomData a;
                a.x = x; a.y = y; a.z = z;
                a.vx = a.vy = a.vz = 0.0;
                a.ax = a.ay = a.az = 0.0;
                atoms.push_back(a);
            }
        }
    }

    // 打印结果（示意）
    cout << "Kept atoms: " << atoms.size() << "\n";
    for (size_t i = 0; i < atoms.size() && i < 10; ++i) {
        cout << i << ": "
             << atoms[i].x << " "
             << atoms[i].y << " "
             << atoms[i].z << "\n";
    }
}
