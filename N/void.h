// #include <string>
// #include <vector>

#include "3.h"
using namespace std;

void readF(const string &configFile, FileName &filename);
void read0(const FileName &filename, Data &data);
void read1(const FileName &filename, parameter1 &p1);
//void read2(const FileName &filename, parameter2 &p2);
void read2(const FileName &filename, parameter2 &p2_WW, parameter2 &p2_BB, parameter2 &p2_WB);
void Dimensionalless(Data &data, parameter1 &pr1, parameter2 &pr2_WW, parameter2 &pr2_BB, parameter2 &pr2_WB);
void RandomV0(Data &data, parameter1 &p1);
//void LJ_potential(Data &data, const parameter1 &pr1, vector<double> &U_atom);
//void BeW_potential(const parameter1 &pr1, const parameter2 &pr2, Data &data, vector<double> &U_atom);
/*void BeW_potential(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
                   */
void BeW_potential2(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
void LJ_evolution(const parameter1 &pr1, Data &Data0,vector<double> &U_atom);

/*void BeW_evolution(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
                   */
void BeW_evolution1(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
void BeW_evolution2(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
void energy(Data &data, const parameter1 &pr1, vector<double> &U_atom);
void InitEnergyFile(const FileName &filename);
void OutputData(const Data &data, const FileName &filename, const parameter1 &pr1);
void OutputEnergy(const Data &data, const FileName &filename, const parameter1 &pr1); 