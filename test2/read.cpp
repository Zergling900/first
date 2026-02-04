#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
//#include <eigen3/Eigen/Dense>

#include "2.h"
// #include "parameter.p1"
//#include "parameter.p2"

//using namespace Eigen;
using namespace std;

void read1(BasicData &data, vector<Atom> &atoms)
{
    data.data_file_name = "error.md3";
    data.T = 0.0;
    data.E = 0.0;
    data.Box_Ln = 0.0;
    data.Box_Lx = 0.0;
    data.Box_Ly = 0.0;
    data.Box_Lz = 0.0;
    data._cell_length_a = 1.0;
    data._cell_length_b = 1.0;
    data._cell_length_c = 1.0;
    data._cell_angle_alpha = 90.0;
    data._cell_angle_beta = 90.0;
    data._cell_angle_gamma = 90.0;
    data.randommultiplier = 0.0;
    data.use_region = 0;
    data.region_ix_min = 0;
    data.region_ix_max = 0;
    data.region_iy_min = 0;
    data.region_iy_max = 0;
    data.region_iz_min = 0;
    data.region_iz_max = 0;
    data.region_in_box = 0;
    data.use_orientation = 0;
    data.orient_h = 0;
    data.orient_k = 0;
    data.orient_l = 1;
    data.orientation_keep_box = 0;
    data.auto_cover_box = 0;
    data.shape_type = "box";
    data.bump_base_ratio = 0.0;
    data.bump_rx = 0.5;
    data.bump_ry = 0.5;
    data.bump_mode = 0.0;
    ifstream fin("parameter.p1");
    if (!fin) {
        cerr << "con't open the file: " << "parameter.p1" << endl;
        return;
    }

    string line, key;
    double value;

    while (getline(fin, line)) {
        if (line.empty() || 
            line[0] == '#' || 
            (line.size() > 1 && line[0] == '/' && line[1] == '/'))
            continue; 

        stringstream ss(line);
        ss >> key;

        if (key == "data_file_name") //string(key) otput file name
        {
            ss >> data.data_file_name;
        }
        else if (key == "shape_type")
        {
            ss >> data.shape_type;
        }
        else
        {
            ss >> value;
            
            if (key == "time") data.T = value;
//            else if (key == "n") data.n =value;
            else if (key == "energy") data.E = value;
            else if (key == "_cell_length_a") data._cell_length_a = value;
            else if (key == "_cell_length_b") data._cell_length_b = value;
            else if (key == "_cell_length_c") data._cell_length_c = value;
            else if (key == "_cell_angle_alpha") data._cell_angle_alpha = value;
            else if (key == "_cell_angle_beta") data._cell_angle_beta = value;
            else if (key == "_cell_angle_gamma") data._cell_angle_gamma = value;
            else if (key == "Box_Ln") data.Box_Ln = value;
            else if (key == "Box_Lx") data.Box_Lx = value;
            else if (key == "Box_Ly") data.Box_Ly = value;
            else if (key == "Box_Lz") data.Box_Lz = value;
            else if (key == "use_region") data.use_region = static_cast<int>(value);
            else if (key == "region_ix_min") data.region_ix_min = static_cast<int>(value);
            else if (key == "region_ix_max") data.region_ix_max = static_cast<int>(value);
            else if (key == "region_iy_min") data.region_iy_min = static_cast<int>(value);
            else if (key == "region_iy_max") data.region_iy_max = static_cast<int>(value);
            else if (key == "region_iz_min") data.region_iz_min = static_cast<int>(value);
            else if (key == "region_iz_max") data.region_iz_max = static_cast<int>(value);
            else if (key == "region_in_box") data.region_in_box = static_cast<int>(value);
            else if (key == "use_orientation") data.use_orientation = static_cast<int>(value);
            else if (key == "orient_h") data.orient_h = static_cast<int>(value);
            else if (key == "orient_k") data.orient_k = static_cast<int>(value);
            else if (key == "orient_l") data.orient_l = static_cast<int>(value);
            else if (key == "orientation_keep_box") data.orientation_keep_box = static_cast<int>(value);
            else if (key == "auto_cover_box") data.auto_cover_box = static_cast<int>(value);
            else if (key == "bump_base_ratio") data.bump_base_ratio = value;
            else if (key == "bump_rx") data.bump_rx = value;
            else if (key == "bump_ry") data.bump_ry = value;
            else if (key == "bump_mode") data.bump_mode = value;
            else if (key == "coordinates") 
            {
                while (getline(fin, line)) 
                {
                    if (line.empty() || 
                        line[0] == '#' || 
                        (line.size() > 1 && line[0] == '/' && line[1] == '/'))
                        continue; 
                    else
                    {
                        Atom atom;
                        stringstream ss(line);
                        ss >> atom.type >> atom.label >> atom.symmetry_multiplicity 
                        >> atom.x0 >> atom.y0 >> atom.z0 >> atom.occupancy;
                        atoms.push_back(atom);
                    }
                }
            }
        }
    }

    fin.close();
}

void read2(Data &data)
{
    ifstream fin("parameter.p2");
    if (!fin) {
        cerr << "con't open the file: " << "parameter.p2" << endl;
        return;
    }

    string line, key;
    double value;

    while (getline(fin, line)) {
        if (line.empty() || 
            line[0] == '#' || 
            (line.size() > 1 && line[0] == '/' && line[1] == '/'))
            continue; 

        stringstream ss(line);
        ss >> key >> value;

        if (key == "x0") data.x = value;
        else if (key == "y0") data.y = value;
        else if (key == "z0") data.z = value;
        else if (key == "xv0") data.vx = value;
        else if (key == "yv0") data.vy = value;
        else if (key == "zv0") data.vz = value;
        else if (key == "dvx0") data.dvx = value;
        else if (key == "dvy") data.dvy = value;
        else if (key == "dvz") data.dvz = value;
    }

    fin.close();
}
