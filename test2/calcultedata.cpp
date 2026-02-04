#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cctype>
//#include <eigen3/Eigen/Dense>

#include "2.h"

//using namespace Eigen;
using namespace std;

//----------------------------------------
// make first unit data
void build(const vector<Atom> &atoms,
           const Data &data1,
           vector<Data> &datas)
{
    datas.clear();
    datas.reserve(atoms.size());

    for (const auto &atom : atoms)
    {
        Data d = data1;
        d.name = atom.type;
        d.x = atom.x0;
        d.y = atom.y0;
        d.z = atom.z0;
        datas.push_back(d);
    }
}

static int box_count(double v, double fallback)
{
    int n = (v > 0.0) ? static_cast<int>(v) : static_cast<int>(fallback);
    if (n < 0)
        n = 0;
    return n;
}

struct Bounds
{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
};

static string to_lower(string s)
{
    for (char &c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

static double clamp01(double v)
{
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

static double clamp_range(double v, double lo, double hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static double det33(const Matrix33 &m)
{
    return
        m.a00 * (m.a11 * m.a22 - m.a12 * m.a21) -
        m.a01 * (m.a10 * m.a22 - m.a12 * m.a20) +
        m.a02 * (m.a10 * m.a21 - m.a11 * m.a20);
}

static bool invert33(const Matrix33 &m, Matrix33 &inv)
{
    const double det = det33(m);
    const double eps = 1e-14;
    if (std::fabs(det) <= eps)
        return false;

    const double inv_det = 1.0 / det;
    inv = Matrix33(
        (m.a11 * m.a22 - m.a12 * m.a21) * inv_det,
        (m.a02 * m.a21 - m.a01 * m.a22) * inv_det,
        (m.a01 * m.a12 - m.a02 * m.a11) * inv_det,

        (m.a12 * m.a20 - m.a10 * m.a22) * inv_det,
        (m.a00 * m.a22 - m.a02 * m.a20) * inv_det,
        (m.a02 * m.a10 - m.a00 * m.a12) * inv_det,

        (m.a10 * m.a21 - m.a11 * m.a20) * inv_det,
        (m.a01 * m.a20 - m.a00 * m.a21) * inv_det,
        (m.a00 * m.a11 - m.a01 * m.a10) * inv_det
    );
    return true;
}

static Bounds bounds_from_parallelepiped(const Matrix31 &o,
                                         const Matrix31 &a,
                                         const Matrix31 &b,
                                         const Matrix31 &c)
{
    Matrix31 p0 = o;
    Matrix31 p1 = o + a;
    Matrix31 p2 = o + b;
    Matrix31 p3 = o + c;
    Matrix31 p4 = o + a + b;
    Matrix31 p5 = o + a + c;
    Matrix31 p6 = o + b + c;
    Matrix31 p7 = o + a + b + c;

    Bounds bd;
    bd.x_min = bd.x_max = p0.a00;
    bd.y_min = bd.y_max = p0.a10;
    bd.z_min = bd.z_max = p0.a20;

    Matrix31 pts[7] = {p1, p2, p3, p4, p5, p6, p7};
    for (const auto &p : pts)
    {
        bd.x_min = std::min(bd.x_min, p.a00);
        bd.x_max = std::max(bd.x_max, p.a00);
        bd.y_min = std::min(bd.y_min, p.a10);
        bd.y_max = std::max(bd.y_max, p.a10);
        bd.z_min = std::min(bd.z_min, p.a20);
        bd.z_max = std::max(bd.z_max, p.a20);
    }
    return bd;
}

static Bounds bounds_from_box(const Matrix33 &Box)
{
    Matrix31 o(0.0, 0.0, 0.0);
    Matrix31 a(Box.a00, Box.a10, Box.a20);
    Matrix31 b(Box.a01, Box.a11, Box.a21);
    Matrix31 c(Box.a02, Box.a12, Box.a22);
    return bounds_from_parallelepiped(o, a, b, c);
}

static bool inside_box(const Matrix31 &r, const Matrix33 &BoxInv)
{
    const double eps = 1e-10;
    Matrix31 f = BoxInv * r;
    return (f.a00 >= -eps && f.a00 <= 1.0 + eps) &&
           (f.a10 >= -eps && f.a10 <= 1.0 + eps) &&
           (f.a20 >= -eps && f.a20 <= 1.0 + eps);
}

static void auto_cover_range_from_box(const Matrix33 &Box,
                                      const Matrix33 &CellInv,
                                      int &ix_min, int &ix_max,
                                      int &iy_min, int &iy_max,
                                      int &iz_min, int &iz_max)
{
    Matrix31 o(0.0, 0.0, 0.0);
    Matrix31 a(Box.a00, Box.a10, Box.a20);
    Matrix31 b(Box.a01, Box.a11, Box.a21);
    Matrix31 c(Box.a02, Box.a12, Box.a22);

    Matrix31 corners[8] = {
        o,
        a,
        b,
        c,
        a + b,
        a + c,
        b + c,
        a + b + c
    };

    Matrix31 f0 = CellInv * corners[0];
    double fx_min = f0.a00, fx_max = f0.a00;
    double fy_min = f0.a10, fy_max = f0.a10;
    double fz_min = f0.a20, fz_max = f0.a20;

    for (int i = 1; i < 8; ++i)
    {
        Matrix31 f = CellInv * corners[i];
        fx_min = std::min(fx_min, f.a00);
        fx_max = std::max(fx_max, f.a00);
        fy_min = std::min(fy_min, f.a10);
        fy_max = std::max(fy_max, f.a10);
        fz_min = std::min(fz_min, f.a20);
        fz_max = std::max(fz_max, f.a20);
    }

    const int pad = 1;
    ix_min = static_cast<int>(std::floor(fx_min)) - pad;
    ix_max = static_cast<int>(std::ceil(fx_max)) + pad;
    iy_min = static_cast<int>(std::floor(fy_min)) - pad;
    iy_max = static_cast<int>(std::ceil(fy_max)) + pad;
    iz_min = static_cast<int>(std::floor(fz_min)) - pad;
    iz_max = static_cast<int>(std::ceil(fz_max)) + pad;
}

static bool inside_shape(const Data &d, const BasicData &data, const Bounds &b)
{
    const double eps = 1e-10;
    const string shape = to_lower(data.shape_type);

    if (shape.empty() || shape == "box")
        return true;

    const double x_c = 0.5 * (b.x_min + b.x_max);
    const double y_c = 0.5 * (b.y_min + b.y_max);
    const double z_min = b.z_min;
    const double z_max = b.z_max;
    const double H = z_max - z_min;
    const double half_x = 0.5 * (b.x_max - b.x_min);
    const double half_y = 0.5 * (b.y_max - b.y_min);
    if (H <= 0.0)
        return true;

    if (shape == "pyramid")
    {
        if (d.z < z_min - eps || d.z > z_max + eps)
            return false;
        const double t = (d.z - z_min) / H;
        const double lim_x = (1.0 - t) * half_x;
        const double lim_y = (1.0 - t) * half_y;
        return (std::fabs(d.x - x_c) <= lim_x + eps) &&
               (std::fabs(d.y - y_c) <= lim_y + eps);
    }

    if (shape == "cone")
    {
        if (d.z < z_min - eps || d.z > z_max + eps)
            return false;
        const double t = (d.z - z_min) / H;
        const double R = std::min(half_x, half_y);
        const double lim_r = (1.0 - t) * R;
        const double dx = d.x - x_c;
        const double dy = d.y - y_c;
        return (dx * dx + dy * dy) <= (lim_r + eps) * (lim_r + eps);
    }

    if (shape == "sphere")
    {
        const double R = 0.5 * H;
        const double z_c = z_min + R;
        const double dx = d.x - x_c;
        const double dy = d.y - y_c;
        const double dz = d.z - z_c;
        return (dx * dx + dy * dy + dz * dz) <= (R + eps) * (R + eps);
    }

    if (shape == "cylinder" || shape == "disk")
    {
        if (d.z < z_min - eps || d.z > z_max + eps)
            return false;
        const double R = std::min(half_x, half_y);
        const double dx = d.x - x_c;
        const double dy = d.y - y_c;
        return (dx * dx + dy * dy) <= (R + eps) * (R + eps);
    }

    if (shape == "bump")
    {
        const double base_ratio = clamp01(data.bump_base_ratio);
        const double base_z = z_min + base_ratio * H;
        const double rx = clamp01(data.bump_rx);
        const double ry = clamp01(data.bump_ry);
        const double half_rx = 0.5 * rx * (b.x_max - b.x_min);
        const double half_ry = 0.5 * ry * (b.y_max - b.y_min);
        const bool in_region = (std::fabs(d.x - x_c) <= half_rx + eps) &&
                               (std::fabs(d.y - y_c) <= half_ry + eps);

        double mode = clamp_range(data.bump_mode, -1.0, 1.0);
        double z_target = base_z;
        if (mode >= 0.0)
            z_target = base_z + mode * (z_max - base_z);
        else
            z_target = base_z + mode * (base_z - z_min);

        if (d.z < z_min - eps)
            return false;
        const double z_limit = in_region ? z_target : base_z;
        return d.z <= z_limit + eps;
    }

    return true;
}

static double vec_norm(const Matrix31 &v)
{
    return std::sqrt(v.a00 * v.a00 + v.a10 * v.a10 + v.a20 * v.a20);
}

static Matrix31 vec_normalize(const Matrix31 &v)
{
    const double n = vec_norm(v);
    if (n <= 0.0)
        return Matrix31(0.0, 0.0, 1.0);
    return Matrix31(v.a00 / n, v.a10 / n, v.a20 / n);
}

static Matrix31 vec_cross(const Matrix31 &a, const Matrix31 &b)
{
    return Matrix31(
        a.a10 * b.a20 - a.a20 * b.a10,
        a.a20 * b.a00 - a.a00 * b.a20,
        a.a00 * b.a10 - a.a10 * b.a00
    );
}

static Matrix33 rotation_align_to_z(const Matrix31 &v)
{
    Matrix31 z = vec_normalize(v);
    Matrix31 up(0.0, 0.0, 1.0);
    if (std::fabs(z.a20) > 0.999)
        up = Matrix31(0.0, 1.0, 0.0);

    Matrix31 x = vec_cross(up, z);
    double xn = vec_norm(x);
    if (xn <= 0.0)
    {
        x = Matrix31(1.0, 0.0, 0.0);
        xn = 1.0;
    }
    x = Matrix31(x.a00 / xn, x.a10 / xn, x.a20 / xn);
    Matrix31 y = vec_cross(z, x);

    return Matrix33(
        x.a00, x.a10, x.a20,
        y.a00, y.a10, y.a20,
        z.a00, z.a10, z.a20
    );
}

Matrix33 ApplyOrientation(const BasicData &data, const Matrix33 &Cell)
{
    if (data.use_orientation == 0)
        return Cell;

    const int h = data.orient_h;
    const int k = data.orient_k;
    const int l = data.orient_l;
    if (h == 0 && k == 0 && l == 0)
        return Cell;

    Matrix31 v_frac(static_cast<double>(h),
                    static_cast<double>(k),
                    static_cast<double>(l));
    Matrix31 v = Cell * v_frac;
    if (vec_norm(v) <= 0.0)
        return Cell;

    Matrix33 R = rotation_align_to_z(v);
    return R * Cell;
}
void calculate(int &n,
               const BasicData &data,
               const Matrix33 &Cell,
               const Matrix33 &BoxCell,
               const vector<Data> &datas,
               vector<Data> &output)
{
    // --------------------------------------------
    output.clear();

    const int nx = box_count(data.Box_Lx, data.Box_Ln);
    const int ny = box_count(data.Box_Ly, data.Box_Ln);
    const int nz = box_count(data.Box_Lz, data.Box_Ln);

    Matrix33 scale(
        nx + 1, 0.0, 0.0,
        0.0, ny + 1, 0.0,
        0.0, 0.0, nz + 1
    );
    Matrix33 Box = BoxCell * scale;
    Bounds b = bounds_from_box(Box);
    Matrix33 BoxInv;
    bool has_box_inv = invert33(Box, BoxInv);
    Matrix33 CellInv;
    bool has_cell_inv = invert33(Cell, CellInv);

    int ix_min = 0, ix_max = nx;
    int iy_min = 0, iy_max = ny;
    int iz_min = 0, iz_max = nz;

    const bool use_auto_cover = (data.auto_cover_box != 0 && has_cell_inv);
    if (use_auto_cover)
    {
        auto_cover_range_from_box(Box, CellInv,
                                  ix_min, ix_max,
                                  iy_min, iy_max,
                                  iz_min, iz_max);
    }
    else if (data.use_region != 0 && data.region_in_box == 0)
    {
        ix_min = data.region_ix_min;
        ix_max = data.region_ix_max;
        iy_min = data.region_iy_min;
        iy_max = data.region_iy_max;
        iz_min = data.region_iz_min;
        iz_max = data.region_iz_max;
    }

    auto clamp_axis = [&](int &min_v, int &max_v, int axis_max) {
        min_v = std::max(0, std::min(min_v, axis_max));
        max_v = std::max(0, std::min(max_v, axis_max));
        if (min_v > max_v)
            std::swap(min_v, max_v);
    };
    if (!use_auto_cover)
    {
        clamp_axis(ix_min, ix_max, nx);
        clamp_axis(iy_min, iy_max, ny);
        clamp_axis(iz_min, iz_max, nz);
    }

    Bounds shape_bounds = bounds_from_box(Box);
    if (data.use_region != 0)
    {
        if (data.region_in_box != 0)
        {
            const double denom_x = (nx + 1.0);
            const double denom_y = (ny + 1.0);
            const double denom_z = (nz + 1.0);

            int rix0 = data.region_ix_min;
            int rix1 = data.region_ix_max;
            int riy0 = data.region_iy_min;
            int riy1 = data.region_iy_max;
            int riz0 = data.region_iz_min;
            int riz1 = data.region_iz_max;
            if (rix0 > rix1) std::swap(rix0, rix1);
            if (riy0 > riy1) std::swap(riy0, riy1);
            if (riz0 > riz1) std::swap(riz0, riz1);

            const double fx0 = clamp_range(rix0 / denom_x, 0.0, 1.0);
            const double fy0 = clamp_range(riy0 / denom_y, 0.0, 1.0);
            const double fz0 = clamp_range(riz0 / denom_z, 0.0, 1.0);
            const double dx = clamp_range((rix1 - rix0 + 1) / denom_x, 0.0, 1.0);
            const double dy = clamp_range((riy1 - riy0 + 1) / denom_y, 0.0, 1.0);
            const double dz = clamp_range((riz1 - riz0 + 1) / denom_z, 0.0, 1.0);

            Matrix31 a(Box.a00, Box.a10, Box.a20);
            Matrix31 b(Box.a01, Box.a11, Box.a21);
            Matrix31 c(Box.a02, Box.a12, Box.a22);
            Matrix31 o = fx0 * a + fy0 * b + fz0 * c;
            Matrix31 ar = dx * a;
            Matrix31 br = dy * b;
            Matrix31 cr = dz * c;
            shape_bounds = bounds_from_parallelepiped(o, ar, br, cr);
        }
        else
        {
            int rix0 = data.region_ix_min;
            int rix1 = data.region_ix_max;
            int riy0 = data.region_iy_min;
            int riy1 = data.region_iy_max;
            int riz0 = data.region_iz_min;
            int riz1 = data.region_iz_max;
            if (rix0 > rix1) std::swap(rix0, rix1);
            if (riy0 > riy1) std::swap(riy0, riy1);
            if (riz0 > riz1) std::swap(riz0, riz1);

            const double dx = static_cast<double>(rix1 - rix0 + 1);
            const double dy = static_cast<double>(riy1 - riy0 + 1);
            const double dz = static_cast<double>(riz1 - riz0 + 1);

            Matrix31 a(Cell.a00, Cell.a10, Cell.a20);
            Matrix31 b(Cell.a01, Cell.a11, Cell.a21);
            Matrix31 c(Cell.a02, Cell.a12, Cell.a22);
            Matrix31 o = (static_cast<double>(rix0) * a) +
                         (static_cast<double>(riy0) * b) +
                         (static_cast<double>(riz0) * c);
            Matrix31 ar = dx * a;
            Matrix31 br = dy * b;
            Matrix31 cr = dz * c;
            shape_bounds = bounds_from_parallelepiped(o, ar, br, cr);
        }
    }

    output.clear();
    for (int ix = ix_min; ix <= ix_max; ++ix)
    {
        for (int iy = iy_min; iy <= iy_max; ++iy)
        {
            for (int iz = iz_min; iz <= iz_max; ++iz)
            {
                if (data.use_region != 0 && data.region_in_box == 0)
                {
                    if (ix < data.region_ix_min || ix > data.region_ix_max ||
                        iy < data.region_iy_min || iy > data.region_iy_max ||
                        iz < data.region_iz_min || iz > data.region_iz_max)
                    {
                        continue;
                    }
                }
                for (const auto &d : datas)
                {
                    Data aaa = d;

                    Matrix31 frac(d.x, d.y, d.z);
                    Matrix31 step(ix, iy, iz);
                    Matrix31 pos = frac + step;
                    //Matrix31 r = pos;
                    Matrix31 r = Cell * pos;
                    //Matrix31 r = Cell_L0 * pos;
                    aaa.x = r.a00;
                    aaa.y = r.a10;
                    aaa.z = r.a20;

                    if (has_box_inv)
                    {
                        if (!inside_box(r, BoxInv))
                            continue;
                    }

                    if (data.use_region != 0 && data.region_in_box != 0 && has_box_inv)
                    {
                        const double eps = 1e-10;
                        Matrix31 f = BoxInv * r;
                        const double vx = f.a00 * (nx + 1.0);
                        const double vy = f.a10 * (ny + 1.0);
                        const double vz = f.a20 * (nz + 1.0);
                        const int ix_cell = static_cast<int>(std::floor(vx + eps));
                        const int iy_cell = static_cast<int>(std::floor(vy + eps));
                        const int iz_cell = static_cast<int>(std::floor(vz + eps));
                        if (ix_cell < data.region_ix_min || ix_cell > data.region_ix_max ||
                            iy_cell < data.region_iy_min || iy_cell > data.region_iy_max ||
                            iz_cell < data.region_iz_min || iz_cell > data.region_iz_max)
                        {
                            continue;
                        }
                    }

                    if (inside_shape(aaa, data, shape_bounds))
                        output.push_back(aaa);
                }
            }
        }
    }
    n = static_cast<int>(output.size());
}

void Output(const int &n, const BasicData &data,
                vector<Data> &output,
                const Matrix33 &BoxCell,const string &filename)
{
    FILE *fp = fopen(filename.c_str(), "w");
    const int nx = box_count(data.Box_Lx, data.Box_Ln);
    const int ny = box_count(data.Box_Ly, data.Box_Ln);
    const int nz = box_count(data.Box_Lz, data.Box_Ln);
    Matrix33 scale(
        nx + 1, 0.0, 0.0,
        0.0, ny + 1, 0.0,
        0.0, 0.0, nz + 1
    );
    Matrix33 Box = BoxCell * scale;
    fprintf(fp, "%d\n", n);
    fprintf(fp, "   time=   %f (fs)  Energy=  %f (eV)\n", data.T, data.E);
    fprintf(fp, "BOX %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
         Box.a00, Box.a10, Box.a20, Box.a01, Box.a11, Box.a21, Box.a02, Box.a12, Box.a22);
    for(const auto &d : output)
    {
        fprintf(fp, " %4s  %18.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
            d.name.c_str(),d.x,d.y,d.z,d.vx,d.vy,d.vz,d.dvx,d.dvy,d.dvz);
    }
    fclose(fp);
}
