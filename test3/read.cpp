#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

#include "3.h"

using namespace std;

//---------------------------------------------------------------------------
void readF(FileName &filename)
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

            if (key == "BasicData_file") filename.BasicData_filename = name;//string(key) otput file name
            else if (key == "Data_file") filename.Data_filename = name;
            else if (key == "Ut_file") filename.Ut_file = name;
            else if (key == "Kt_file") filename.Kt_file = name;
            //else if (key == "other1") filename.other1 = name;
            //else if (key == "other2") filename.other2 = name;
    }
};

void read0(const FirstData &data0,const FileName &filename, Data &data)
{
    ifstream fin(filename.BasicData_filename);
    if (!fin) {
        cerr << "con't open the file: " << filename.BasicData_filename << endl;
        return;
    }

    string line, key;
    double value;

}