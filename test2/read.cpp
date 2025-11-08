#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <eigen3/Eigen/Dense>

#include "2.h"
#include "parameter.p1"
//#include "parameter.p2"

using namespace Eigen;
using namespace std;


#include <iostream>

void read(BasicData &data)
{
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
        ss >> key >> value;

        if (key == "energy") data.E = value;
        else if (key == "time") data.T = value;

        else if (key == "Cell_La") data.Cell_La = value;
        else if (key == "Cell_Lb") data.Cell_Lb = value;
        else if (key == "Cell_Lc") data.Cell_Lc = value;
        else if (key == "cell_angle_alpha") data.cell_angle_alpha = value;
        else if (key == "cell_angle_beta") data.cell_angle_beta = value;
        else if (key == "cell_angle_gamma") data.cell_angle_gamma = value;
        else if (key == "Box_Ln") data.Box_Ln = value;
        else if (key == "Box_Lx") data.Box_Lx = value;
        else if (key == "Box_Ly") data.Box_Ly = value;
        else if (key == "Box_Lz") data.Box_Lz = value;
        else if (key == "coordinates") 
    }

    fin.close();
}