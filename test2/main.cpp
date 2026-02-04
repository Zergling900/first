#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
//#include <eigen3/Eigen/Dense>

#include "2.h"

//using namespace Eigen;
using namespace std;
//----------------------------------------------
//main
//----------------------------------------------

// main.cpp
#include "2.h"

int main()
{
    BasicData data;
    vector<Atom> atoms;

    Data data1;
    vector<Data> datas;
    vector<Data> output;

    cout << "STEP 1: read1\n";
    read1(data, atoms);

    cout << "STEP 2: read2\n";
    read2(data1);

    cout << "STEP 3: build\n";
    build(atoms, data1, datas);

    cout << "STEP 4: CellVector\n";
    Matrix33 Cell_L0 = CellVector_L0(data);
    Matrix33 Cell_A = CellVector_A(data);
    Matrix33 Cell_base = CellVector_C(Cell_L0, Cell_A);
    Matrix33 Cell_oriented = ApplyOrientation(data, Cell_base);
    Matrix33 BoxCell = Cell_oriented;
    if (data.use_orientation == 0)
        BoxCell = Cell_base;
    else if (data.orientation_keep_box != 0)
        BoxCell = Cell_base;

    cout << "STEP 5: calculate\n";
    int n = 0;
    calculate(n, data, Cell_oriented, BoxCell, datas, output);

    cout << "STEP 6: Output\n";
    Output(n, data, output, BoxCell, data.data_file_name);

    cout << "DONE\n";
    return 0;
}
