#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>

//-----------------
struct Atom {
    double x, y, z;
};

struct AtomCmp {
    bool operator()(const Atom& a, const Atom& b) const {
        auto q = [](double v) {
            return std::llround(v * 1e6); // 量化到 1e-6
        };
        if (q(a.x) != q(b.x)) return q(a.x) < q(b.x);
        if (q(a.y) != q(b.y)) return q(a.y) < q(b.y);
        return q(a.z) < q(b.z);
    }
};

int main() {
    //
    const double a0 = 2.0;
    const int nx = 3;
    const int ny = 3;
    const int nz = 3;
    const bool add_random_shift = true;

    // 
    double shift_x = 0.0, shift_y = 0.0, shift_z = 0.0;
    if (add_random_shift) {
        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> uni(0.0, 0.1 * a0); 

        shift_x = uni(rng);
        shift_y = uni(rng);
        shift_z = uni(rng);
    }

    //bcc basic
    const std::vector<std::array<double,3>> basis = 
    {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };

    //
    std::set<Atom, AtomCmp> atoms;

    // 
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {

                for (auto &b : basis) {
                    Atom a;
                    a.x = (ix + b[0]) * a0 + shift_x;
                    a.y = (iy + b[1]) * a0 + shift_y;
                    a.z = (iz + b[2]) * a0 + shift_z;
                    atoms.insert(a); 
                }
            }
        }
    }

    // output
    std::ofstream fout("atoms.xyz");
    fout << atoms.size() << "\n";
    fout << "BCC lattice a0=" << a0 
         << " nx=" << nx << " ny=" << ny << " nz=" << nz << "\n";

    for (const auto &a : atoms) {
        fout << "W "
             << std::fixed << std::setprecision(8)
             << std::setw(12) << a.x << " "
             << std::setw(12) << a.y << " "
             << std::setw(12) << a.z << "\n";
    }

    fout.close();

    std::cout << "Generated " << atoms.size() 
              << " unique atoms. Wrote atoms.xyz\n";
    std::cout << "Random shift = (" 
              << shift_x << ", " << shift_y << ", " << shift_z << ")\n";

    return 0;
}
