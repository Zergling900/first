// #include <string>
// #include <vector>

#include "3.h"
using namespace std;

void readF(FileName &filename);
void read0(const FileName &filename, Data &data);
void read1(const FileName &filename, parameter1 &p1);
void RandomV0(Data &Data0, const parameter1 &p1);
void LJ_potential(Data &data, const parameter1 &pr1, vector<double> &U_atom);
void evolution(Data &data, const parameter1 &pr1);
void energy(Data &data, const parameter1 &pr1, vector<double> &U_atom);
void InitEnergyFile(const FileName &filename);
void Output(const Data &data, const FileName &filename, const parameter1 &p1);