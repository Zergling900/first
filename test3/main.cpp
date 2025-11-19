#include <iostream>

#include "3.h"
#include "void.h"

int main()
{

    FileName filename;
    Data data;
    parameter1 p1;
    vector<double> U_atom;

    //read basic data*******************************************
    cout << "STEP 1: readF\n";
    readF(filename);

    cout << "STEP 2: read0\n";
    read0(filename, data);

    cout << "STEP 3: read1\n";
    read1(filename, p1);

    cout << "STEP 4: InitEnergyFile\n";
    InitEnergyFile(filename);
    //***********************************************************

    //initial***************************************************
    cout << "STEP 5: RandomV0\n";
    RandomV0(data, p1);

    cout << "STEP 6: First_LJ_potential\n";
    LJ_potential(data, p1, U_atom);
    Output(data, filename,p1);
    //***********************************************************

    //evolution*************************************************
    cout << "STEP 6: evolution\n";
    for(int i = 0; i < p1.steps; i++)
    {
        evolution(data, p1);
        energy(data, p1, U_atom);
        Output(data, filename, p1);
        data.T += p1.dt;
    }
    //***********************************************************

    return 0;
}