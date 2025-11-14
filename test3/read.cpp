#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

#include "3.h"

using namespace std;

//---------------------------------------------------------------------------
void readF(FileData &filename)
{
    ifstream fin("FileName.txt");
    if (!fin) {
        cerr << "con't open the file: " << "FileName.txt" << endl;
        return;
    }

    string line, key, name;

    while (getline(fin, line)) {
        if (line.empty() || 
            line[0] == '#' || 
            (line.size() > 1 && line[0] == '/' && line[1] == '/'))
            continue; 

        stringstream ss(line);
        ss >> key >> name;

            if (key == "BasicData_file") filename.BasicData_file = name;//string(key) otput file name
            else if (key == "Data_file") filename.Data_file = name;
            else if (key == "Ut_file") filename.Ut_file = name;
            else if (key == "Kt_file") filename.Kt_file = name;
    }
};
