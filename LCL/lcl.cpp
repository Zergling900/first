#include <iostream>
#include <random>
#include <cmath>
#include <vector>

#include "3.h"

void lcl(const Data&data,const parameter1 &pr1, const parameter2 &pr2)
{
    Matrix31 Box ={data.Box.a00, data.Box.a11, data.Box.a22};
    double cutoff = pr2.R;
    int atom_in_cell[M][S];
};