#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

#include "3.h"

//33 * 33
Matrix33 operator*(const Matrix33 &A,const Matrix33 &B)
{
    return Matrix33(
        A.a00 * B.a00 + A.a01 * B.a10 + A.a02 * B.a20,    A.a00 * B.a01 + A.a01 * B.a11 + A.a02 * B.a21,    A.a00 * B.a02 + A.a01 * B.a12 + A.a02 * B.a22,
        A.a10 * B.a00 + A.a11 * B.a10 + A.a12 * B.a20,    A.a10 * B.a01 + A.a11 * B.a11 + A.a12 * B.a21,    A.a10 * B.a02 + A.a11 * B.a12 + A.a12 * B.a22,
        A.a20 * B.a00 + A.a21 * B.a10 + A.a22 * B.a20,    A.a20 * B.a01 + A.a21 * B.a11 + A.a22 * B.a21,    A.a20 * B.a02 + A.a21 * B.a12 + A.a22 * B.a22
    );
}
//---------------------------------------------------------------------------------------------------------------
//33 +33
Matrix33 operator+(const Matrix33 &A,const Matrix33 &B)
{
    return Matrix33(
        A.a00 + B.a00, A.a01 + B.a01, A.a02 + B.a02,
        A.a10 + B.a10, A.a11 + B.a11, A.a12 + B.a12,
        A.a20 + B.a20, A.a21 + B.a21, A.a22 + B.a22
    );
}
//---------------------------------------------------------------------------------------------------------------
//33 -33
Matrix33 operator-(const Matrix33 &A,const Matrix33 &B)
{
    return Matrix33(
        A.a00 - B.a00, A.a01 - B.a01, A.a02 - B.a02,
        A.a10 - B.a10, A.a11 - B.a11, A.a12 - B.a12,
        A.a20 - B.a20, A.a21 - B.a21, A.a22 - B.a22
    );
}
//---------------------------------------------------------------------------------------------------------------
//33 *31
Matrix31 operator*(const Matrix33 &A,const Matrix31 &B)
{
    return Matrix31(
        A.a00 * B.a00 + A.a01 * B.a10 + A.a02 * B.a20,
        A.a10 * B.a00 + A.a11 * B.a10 + A.a12 * B.a20,
        A.a20 * B.a00 + A.a21 * B.a10 + A.a22 * B.a20
    );
}
//---------------------------------------------------------------------------------------------------------------
//31 +31
Matrix31 operator+(const Matrix31 &A,const Matrix31 &B)
{
    return Matrix31(
        A.a00 + B.a00, A.a10 + B.a10, A.a20 + B.a20
    );
}
//---------------------------------------------------------------------------------------------------------------
//31 - 31
Matrix31 operator-(const Matrix31 &A,const Matrix31 &B)
{
    return Matrix31(
        A.a00 - B.a00, A.a10 - B.a10, A.a20 - B.a20
    );
}
//---------------------------------------------------------------------------------------------------------------
//k*31
Matrix31 operator*(double k, const Matrix31 &v)
{
    return Matrix31(k * v.a00, k * v.a10, k * v.a20);
}
//---------------------------------------------------------------------------------------------------------------
//31*k
Matrix31 operator*(const Matrix31 &v, double k)
{
    return Matrix31(k * v.a00, k * v.a10, k * v.a20);
}
//---------------------------------------------------------------------------------------------------------------
//k*33
Matrix33 operator*(double k, const Matrix33 &M)
{
    return Matrix33(
        k * M.a00, k * M.a01, k * M.a02,
        k * M.a10, k * M.a11, k * M.a12,
        k * M.a20, k * M.a21, k * M.a22
    );
}
//---------------------------------------------------------------------------------------------------------------
//33*k
Matrix33 operator*(const Matrix33 &M, double k)
{
    return k * M; 
}
