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
// BasicData_filename,Data_filename,Ut_file,Kt_file
void readF(const std::string &configFile, FileName &fn)
{
    std::ifstream fin(configFile);
    if (!fin)
    {
        std::cerr << "can't open file: " << configFile << std::endl;
        return;
    }

    string key, x, line;

    while (getline(fin, line))
    {
        if (line.empty() ||
            line[0] == '#' ||
            (line.size() > 1 && line[0] == '/' && line[1] == '/'))
            continue;

        std::stringstream ss(line);
        ss >> key >> x;
        // parameter                           parameter.p1
        // BCC

        // First_Data_file                          /BasicData/BCC.md3
        // Data_file                                /Data/BCC.md3
        // Et_file                             /Data/BCC.Et.md3
        if (key == "")
            continue; // skip space
        if (x == "")
            continue; // 

        if (key == "parameter")
        {
            fn.parameter_filename = x;
        }
        else if (key == "parameter2")
        {
            fn.parameter2_filename = x;
        }
        else if (key == "First_Data_file")
        {
            fn.BasicData_filename = x;
        }
        else if (key == "Data_file")
        {
            fn.Data_filename = x;
        }
        else if (key == "Et_file")
        {
            fn.Et_file = x;
        }
    }
}

void read0(const FileName &filename, Data &data)
{
    ifstream fin(filename.BasicData_filename);
    if (!fin)
    {
        cerr << "Can't open file: " << filename.BasicData_filename << endl;
        return;
    }

    // 1. read n
    fin >> data.n;
    string dummy;
    getline(fin, dummy);

    // 2. read time & energy
    string tmp;
    fin >> tmp >> data.T >> tmp >> tmp >> data.E;
    getline(fin, dummy);

    // 3. read Box
    fin >> tmp;
    fin >> data.Box.a00 >> data.Box.a01 >> data.Box.a02 >> data.Box.a10 >> data.Box.a11 >> data.Box.a12 >> data.Box.a20 >> data.Box.a21 >> data.Box.a22;

    // 4. read atoms
    data.atoms.resize(data.n);

    for (int i = 0; i < data.n; ++i)
    {
        Atom &a = data.atoms[i];
        fin >> a.name >> a.r.a00 >> a.r.a10 >> a.r.a20 >> a.p.a00 >> a.p.a10 >> a.p.a20 >> a.f.a00 >> a.f.a10 >> a.f.a20;
    }
}

void read1(const FileName &filename, parameter1 &p1)
{
    ifstream fin(filename.parameter_filename);
    if (!fin)
    {
        std::cerr << "can't open file: " << filename.parameter_filename << std::endl;
        return;
    }

    string key, line;
    double value;
    int i;
    while (getline(fin, line))
    {
        stringstream ss(line);

        ss >> key >> value;

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
        else if (key == "T")
        {
            p1.T = value;
        }
        else if (key == "kb")
        {
            p1.kb = value;
        }
        else if (key == "epsilon")
        {
            p1.epsilon = value;
        }
        else if (key == "sigma")
        {
            p1.sigma = value;
        }
        else if (key == "mass")
        {
            p1.m = value;
        }
    }

    fin.close();
};

void read2(const FileName &filename, parameter2 &pr2)
{
    ifstream fin(filename.parameter2_filename);
    if (!fin)
    {
        std::cerr << "can't open file: " << filename.parameter2_filename << std::endl;
        return;
    }

    string key, line;
    double value;
    bool reading_abop = false;

    while (getline(fin, line))
    {
        // ignore //
        size_t pos = line.find("//");
        if (pos != string::npos)
            line = line.substr(0, pos);

        stringstream ss(line);

        if (!(ss >> key))
            continue;

        // ----------------------------------------------------------------

        if (key == "ABOP")
        {
            ss >> key; // Abb / Aww / Abw
            cout << "Using parameter set: " << key << endl;
            reading_abop = true;
        }
        else if (reading_abop)
        {
            if (key == "D0")
                ss >> pr2.D0;
            else if (key == "r0")
                ss >> pr2.r0;
            else if (key == "beta")
                ss >> pr2.beta;
            else if (key == "S")
                ss >> pr2.S;
            else if (key == "gamma")
                ss >> pr2.gamma;
            else if (key == "c")
                ss >> pr2.c;
            else if (key == "d")
                ss >> pr2.d;
            else if (key == "h")
                ss >> pr2.h;
            else if (key == "R")
                ss >> pr2.R;
            else if (key == "D")
                ss >> pr2.D;
            else if (key == "2mu")
                ss >> pr2.mu;
            else if (key == "rf")
                ss >> pr2.rf;
            else if (key == "bf")
                ss >> pr2.bf;
        }
    }

    fin.close();
    // pr2.mu = 0.5 * pr2.mu;
};