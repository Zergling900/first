#pragma once

#include <iostream>

void read(const parameter &parameter.p1, BasicData &data)
{
    ifstream fin(filename);
    if (!fin) {
        cerr << "无法打开文件: " << filename << endl;
        return;
    }

    string line, key;
    double value;

    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#' || (line.size() > 1 && line[0] == '/' && line[1] == '/'))
            continue; // 跳过注释行

        stringstream ss(line);
        ss >> key >> value;

        if (key == "energy") data.energy = value;
        else if (key == "Box_Lx") data.Box_Lx = value;
        else if (key == "Box_Ly") data.Box_Ly = value;
        else if (key == "Box_Lz") data.Box_Lz = value;
        else if (key == "a0") data.a0 = value;
    }

    fin.close();
}