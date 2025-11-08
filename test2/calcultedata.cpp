#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <eigen3/Eigen/Dense>


void calculate(FirstMolecularData &MolecularData, BasicData &data, BCCData &Data)
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