#include <iostream>
// #include <random>
#include <sstream>
#include <string>
// #include <cmath>
// #include <fstream>
#include <vector>
// #include <iomanip>

#include "3.h"

using namespace std;

//---------------------------------------------------------------------------
void read0(const FileName &filename, Data &data, std::vector<Data> &atoms)
{
    ifstream fin(filename.BasicData_filename);
    if (!fin)
    {
        std::cerr << "can't open file: " << filename.BasicData_filename << std::endl;
        return;
    }

    // 1. n
    fin >> data.n;
    string dummy;
    getline(fin, dummy); // skip line

    // 2. time & Energy
    string tmp;
    fin >> tmp;          // "time="
    fin >> data.T;       // T
    fin >> tmp;          // "(fs)"
    fin >> tmp;          // "Energy="
    fin >> data.E;       // E
    getline(fin, dummy); // skip line

    // 3. BOX
    fin >> tmp; // "BOX"
    fin >> data.Box.a00 >> data.Box.a10 >> data.Box.a20 >> data.Box.a01 >> data.Box.a11 >> data.Box.a21 >> data.Box.a02 >> data.Box.a12 >> data.Box.a22;

    // 4. atoms
    atoms.clear();
    atoms.reserve(data.n);

    for (int i = 0; i < data.n - 1; ++i)
    {
        Data d;
        fin >> d.name >> d.x >> d.y >> d.z >> d.vx >> d.vy >> d.vz >> d.dvx >> d.dvy >> d.dvz;
        atoms.push_back(d);
    }
}

void read1(parameter1 &p1)
{
    ifstream fin("parameter.p1");
    if (!fin)
    {
        std::cerr << "can't open file: " << "parameter.p1" << std::endl;
        return;
    }

    string key;
    double value;
    int i;
    while (fin >> key >> value)
    {
        if (key == "dt")
        {
            p1.dt = value;
        }
        else if (key == "steps")
        {
            p1.steps = static_cast<int>(value);
        }
        else if (key == "steps_space")
        {
            p1.steps_space = static_cast<int>(value);
        }
    }

    fin.close();
}