#include <iostream>

#include "3.h"
#include "void.h"

int main()
{

    FileName filename;
    Data data;
    parameter1 pr1;
    parameter2 pr2;
    vector<double> U_atom;

    // read basic data*******************************************
    cout << "STEP 1: readF\n";

    readF("FileName.txt", filename);

    cout << "STEP 2: read0\n";
    read0(filename, data);

    cout << "STEP 3: read1\n";
    read1(filename, pr1);

    cout << "STEP 4: read2\n";
    read2(filename, pr2);

    cout << "STEP 4: InitEnergyFile\n";
    InitEnergyFile(filename);
    //***********************************************************

    // initial***************************************************
    cout << "STEP 6: RandomV0\n";
    RandomV0(data, pr1);

    cout << "STEP 7: First_potential\n";
    //LJ_potential(data, pr1, U_atom);
    BeW_potential(pr1, pr2, data, U_atom);
    energy(data, pr1, U_atom);
    
    OutputData(data, filename, pr1);
    OutputEnergy(data, filename, pr1);
    //OutputData(data, filename,pr1);
    //***********************************************************

    // evolution*************************************************
    cout << "STEP 8: evolution\n";
    for (int i = 0; i < pr1.steps; i++)
    {
        data.T += pr1.dt;
        printf("calculating time = %f\n", data.T);
        //LJ_evolution(pr1, data, U_atom);
        BeW_evolution(pr1, pr2, data, U_atom);
        // energy(data, pr1, U_atom); in evolution

        if ((i+1) % pr1.steps_space == 0)
        {
            OutputData(data, filename, pr1);
            OutputEnergy(data, filename, pr1);
        }
        //data.T += pr1.dt;
    }
    //***********************************************************

    return 0;
}