#include <iostream>
#include <cmath>
#include <vector>

#include "3.h"

void Dimensionalless(Data &data, parameter1 &pr1, parameter2 &pr2_WW, parameter2 &pr2_BB, parameter2 &pr2_WB)
{
    double E_0 = pr1.E0;
    pr1.kb = pr1.kb / E_0;
    
    pr2_WW.D0 = pr2_WW.D0 / pr1.E0;
    pr2_BB.D0 = pr2_BB.D0 / pr1.E0;
    pr2_WB.D0 = pr2_WB.D0 / pr1.E0;
}