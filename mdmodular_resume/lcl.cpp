#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "3.h"
#include "void.h"

namespace
{
inline int WrapIndex(int idx, int m)
{
    if (m <= 0)
        return 0;
    idx %= m;
    if (idx < 0)
        idx += m;
    return idx;
}

inline double WrapDelta(double d, double L, double Lh)
{
    if (d < -Lh)
        d += L;
    else if (d >= Lh)
        d -= L;
    return d;
}

bool NeedVerletRebuild(const Data &data, const Cell_List &cl)
{
    const int n = data.n;
    if (!cl.verlet_valid)
        return true;
    if (cl.verlet_ref_n != n)
        return true;
    if (static_cast<int>(cl.verlet_ref_pos.size()) != n)
        return true;
    if (cl.verlet_rebuild_limit <= 0.0)
        return true;

    const double Lx = data.Box.a00;
    const double Ly = data.Box.a11;
    const double Lz = data.Box.a22;
    const double Lxh = 0.5 * Lx;
    const double Lyh = 0.5 * Ly;
    const double Lzh = 0.5 * Lz;
    const double limit2 = cl.verlet_rebuild_limit * cl.verlet_rebuild_limit;

    for (int i = 0; i < n; ++i)
    {
        const Matrix31 &r = data.atoms[i].r;
        const Matrix31 &rr = cl.verlet_ref_pos[i];
        const double dx = WrapDelta(r.a00 - rr.a00, Lx, Lxh);
        const double dy = WrapDelta(r.a10 - rr.a10, Ly, Lyh);
        const double dz = WrapDelta(r.a20 - rr.a20, Lz, Lzh);
        if (dx * dx + dy * dy + dz * dz > limit2)
            return true;
    }
    return false;
}

void BuildCellOccupancy(const Data &data, Cell_List &cl)
{
    const int n = data.n;
    const int Mx = cl.Mx;
    const int My = cl.My;
    const int Mz = cl.Mz;
    const double Wx = cl.Wx;
    const double Wy = cl.Wy;
    const double Wz = cl.Wz;

    if (static_cast<int>(cl.num_in_cell_cache.size()) != cl.cell_num)
        cl.num_in_cell_cache.resize(cl.cell_num);
    if (static_cast<int>(cl.offset_cache.size()) != cl.cell_num)
        cl.offset_cache.resize(cl.cell_num);

    std::fill(cl.num_in_cell_cache.begin(), cl.num_in_cell_cache.end(), 0);
    std::vector<int> &num_in_cell = cl.num_in_cell_cache;

    for (int i = 0; i < n; ++i)
    {
        int mx = static_cast<int>(std::floor(data.atoms[i].r.a00 / Wx));
        int my = static_cast<int>(std::floor(data.atoms[i].r.a10 / Wy));
        int mz = static_cast<int>(std::floor(data.atoms[i].r.a20 / Wz));
        mx = WrapIndex(mx, Mx);
        my = WrapIndex(my, My);
        mz = WrapIndex(mz, Mz);
        const int cid = mx + my * Mx + mz * Mx * My;
        cl.Cell[i] = cid;
        ++num_in_cell[cid];
    }

    cl.cell_offset[0] = 0;
    for (int c = 0; c < cl.cell_num; ++c)
        cl.cell_offset[c + 1] = cl.cell_offset[c] + num_in_cell[c];

    std::vector<int> &offset1 = cl.offset_cache;
    for (int c = 0; c < cl.cell_num; ++c)
        offset1[c] = cl.cell_offset[c];

    for (int i = 0; i < n; ++i)
    {
        const int c = cl.Cell[i];
        const int pos = offset1[c]++;
        cl.atom_indices[pos] = i;
    }
}

void BuildCellNeighborStencil(Cell_List &cl)
{
    const int cell_num = cl.cell_num;
    const int Mx = cl.Mx;
    const int My = cl.My;
    const int Mz = cl.Mz;

    cl.cell_neighbor_offset.resize(cell_num + 1);
    cl.cell_neighbor_cells.clear();
    cl.cell_neighbor_cells.reserve(static_cast<size_t>(cell_num) * 27u);

    for (int cid = 0; cid < cell_num; ++cid)
    {
        cl.cell_neighbor_offset[cid] = static_cast<int>(cl.cell_neighbor_cells.size());

        const int cz = cid / (Mx * My);
        const int cc = cid - cz * (Mx * My);
        const int cy = cc / Mx;
        const int cx = cc - cy * Mx;

        int neighbor_cells[27];
        int neighbor_count = 0;
        for (int dcz = -1; dcz <= 1; ++dcz)
        {
            for (int dcy = -1; dcy <= 1; ++dcy)
            {
                for (int dcx = -1; dcx <= 1; ++dcx)
                {
                    const int cz2 = WrapIndex(cz + dcz, Mz);
                    const int cy2 = WrapIndex(cy + dcy, My);
                    const int cx2 = WrapIndex(cx + dcx, Mx);
                    const int nc = cx2 + cy2 * Mx + cz2 * Mx * My;
                    bool duplicate = false;
                    for (int t = 0; t < neighbor_count; ++t)
                    {
                        if (neighbor_cells[t] == nc)
                        {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate)
                        neighbor_cells[neighbor_count++] = nc;
                }
            }
        }

        for (int t = 0; t < neighbor_count; ++t)
            cl.cell_neighbor_cells.push_back(neighbor_cells[t]);
    }
    cl.cell_neighbor_offset[cell_num] = static_cast<int>(cl.cell_neighbor_cells.size());
}

void BuildVerletNeighborList(const Data &data, Cell_List &cl)
{
    const int n = data.n;
    const double Lx = data.Box.a00;
    const double Ly = data.Box.a11;
    const double Lz = data.Box.a22;
    const double Lxh = 0.5 * Lx;
    const double Lyh = 0.5 * Ly;
    const double Lzh = 0.5 * Lz;
    const double rlist2 = cl.verlet_cutoff2;

    if (static_cast<int>(cl.verlet_count_cache.size()) != n)
        cl.verlet_count_cache.resize(n);
    if (static_cast<int>(cl.verlet_cursor_cache.size()) != n)
        cl.verlet_cursor_cache.resize(n);
    std::fill(cl.verlet_count_cache.begin(), cl.verlet_count_cache.end(), 0);

    if (static_cast<int>(cl.cell_neighbor_offset.size()) != cl.cell_num + 1)
        BuildCellNeighborStencil(cl);

    std::vector<int> &pair_i_cache = cl.verlet_pair_i_cache;
    std::vector<int> &pair_j_cache = cl.verlet_pair_j_cache;
    pair_i_cache.clear();
    pair_j_cache.clear();

    size_t reserve_pairs = cl.verlet_neighbors.empty() ? 0u : static_cast<size_t>(cl.verlet_neighbors.size() / 2);
    const size_t min_pairs = static_cast<size_t>(n > 0 ? n : 1);
    if (reserve_pairs < min_pairs)
        reserve_pairs = min_pairs;
    pair_i_cache.reserve(reserve_pairs);
    pair_j_cache.reserve(reserve_pairs);

    auto register_pair_if_inside_cutoff = [&](int i, int j)
    {
        const double dx = WrapDelta(data.atoms[j].r.a00 - data.atoms[i].r.a00, Lx, Lxh);
        const double dy = WrapDelta(data.atoms[j].r.a10 - data.atoms[i].r.a10, Ly, Lyh);
        const double dz = WrapDelta(data.atoms[j].r.a20 - data.atoms[i].r.a20, Lz, Lzh);
        const double rij2 = dx * dx + dy * dy + dz * dz;
        if (rij2 < rlist2)
        {
            ++cl.verlet_count_cache[i];
            ++cl.verlet_count_cache[j];
            pair_i_cache.push_back(i);
            pair_j_cache.push_back(j);
        }
    };

    for (int i = 0; i < n; ++i)
    {
        const int c_i = cl.Cell[i];
        const int h_begin = cl.cell_neighbor_offset[c_i];
        const int h_end = cl.cell_neighbor_offset[c_i + 1];
        for (int h = h_begin; h < h_end; ++h)
        {
            const int nc = cl.cell_neighbor_cells[h];
            const int begin = cl.cell_offset[nc];
            const int end = cl.cell_offset[nc + 1];
            for (int jj = begin; jj < end; ++jj)
            {
                const int j = cl.atom_indices[jj];
                if (j <= i)
                    continue;
                register_pair_if_inside_cutoff(i, j);
            }
        }
    }

    cl.verlet_offset.resize(n + 1);
    cl.verlet_offset[0] = 0;
    for (int i = 0; i < n; ++i)
        cl.verlet_offset[i + 1] = cl.verlet_offset[i] + cl.verlet_count_cache[i];

    cl.verlet_neighbors.resize(cl.verlet_offset[n]);
    for (int i = 0; i < n; ++i)
        cl.verlet_cursor_cache[i] = cl.verlet_offset[i];

    const int pair_count = static_cast<int>(pair_i_cache.size());
    for (int idx = 0; idx < pair_count; ++idx)
    {
        const int i = pair_i_cache[idx];
        const int j = pair_j_cache[idx];
        cl.verlet_neighbors[cl.verlet_cursor_cache[i]++] = j;
        cl.verlet_neighbors[cl.verlet_cursor_cache[j]++] = i;
    }

    cl.verlet_ref_pos.resize(n);
    for (int i = 0; i < n; ++i)
        cl.verlet_ref_pos[i] = data.atoms[i].r;
    cl.verlet_ref_n = n;
    cl.verlet_valid = true;
    ++cl.verlet_build_count;
}
} // namespace

void lcl0(Data &data, Cell_List &cl, const parameter1 &pr1,
          const parameter2 &pr2_WW,
          const parameter2 &pr2_BB,
          const parameter2 &pr2_WB)
{
    (void)pr1;
    const int n = data.n;

    const double cutoff_WW = pr2_WW.R + pr2_WW.D;
    const double cutoff_BB = pr2_BB.R + pr2_BB.D;
    const double cutoff_WB = pr2_WB.R + pr2_WB.D;
    const double pair_cutoff = std::max(cutoff_WW, std::max(cutoff_BB, cutoff_WB));
    if (pair_cutoff <= 0.0 || data.Box.a00 <= 0.0 || data.Box.a11 <= 0.0 || data.Box.a22 <= 0.0)
    {
        std::cerr << "[lcl0] Invalid box or cutoff (R+D).\n";
        cl = Cell_List{};
        return;
    }

    if (cl.verlet_skin <= 0.0)
        cl.verlet_skin = 0.35;
    cl.verlet_cutoff = pair_cutoff + cl.verlet_skin;
    cl.verlet_cutoff2 = cl.verlet_cutoff * cl.verlet_cutoff;
    cl.verlet_rebuild_limit = 0.5 * cl.verlet_skin;

    cl.Mx = static_cast<int>(std::floor(data.Box.a00 / cl.verlet_cutoff));
    cl.My = static_cast<int>(std::floor(data.Box.a11 / cl.verlet_cutoff));
    cl.Mz = static_cast<int>(std::floor(data.Box.a22 / cl.verlet_cutoff));
    if (cl.Mx < 1)
        cl.Mx = 1;
    if (cl.My < 1)
        cl.My = 1;
    if (cl.Mz < 1)
        cl.Mz = 1;

    cl.Wx = data.Box.a00 / cl.Mx;
    cl.Wy = data.Box.a11 / cl.My;
    cl.Wz = data.Box.a22 / cl.Mz;

    cl.cell_num = cl.Mx * cl.My * cl.Mz;
    cl.Cell.resize(n);
    cl.cell_offset.resize(cl.cell_num + 1);
    cl.atom_indices.resize(n);
    cl.cell_neighbor_offset.resize(cl.cell_num + 1);
    cl.cell_neighbor_cells.clear();
    cl.num_in_cell_cache.resize(cl.cell_num);
    cl.offset_cache.resize(cl.cell_num);
    cl.verlet_count_cache.resize(n);
    cl.verlet_cursor_cache.resize(n);
    cl.verlet_pair_i_cache.clear();
    cl.verlet_pair_j_cache.clear();
    cl.verlet_ref_pos.resize(n);
    BuildCellNeighborStencil(cl);
    cl.verlet_valid = false;
    cl.verlet_ref_n = n;
}

void lcl1(Data &data, Cell_List &cl, const parameter1 &pr1, const parameter2 &pr2, vector<double> &U_atom)
{
    (void)pr1;
    (void)pr2;
    (void)U_atom;

    if (cl.cell_num <= 0 || cl.Mx <= 0 || cl.My <= 0 || cl.Mz <= 0)
    {
        std::cerr << "[lcl1] Invalid cell list.\n";
        return;
    }

    if (!NeedVerletRebuild(data, cl))
        return;

    BuildCellOccupancy(data, cl);
    BuildVerletNeighborList(data, cl);
}
