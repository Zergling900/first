#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

#include "2.h"
struct Matrix33
{
    double
    a00,a01,a02,
    a10,a11,a12,
    a20,a21,a22;
    Matrix33(double m00, double m01, double m02,
             double m10, double m11, double m12,
             double m20, double m21, double m22)
        : a00(m00), a01(m01), a02(m02),
          a10(m10), a11(m11), a12(m12),
          a20(m20), a21(m21), a22(m22) {}
};

struct Matrix31
{
    double
    x,y,z;
    Matrix31(double m00, double m01, double m02)
        : x(m00), y(m01), z(m02) {}
};


void MatrixTimes3333(const Matrix33 &A,const Matrix33 &B, Matrix33 &C)
{
    C.a00 = A.a00 * B.a00 + A.a01 * B.a10 + A.a02 * B.a20;
    C.a01 = A.a00 * B.a01 + A.a01 * B.a11 + A.a02 * B.a21;
    C.a02 = A.a00 * B.a02 + A.a01 * B.a12 + A.a02 * B.a22;
    C.a10 = A.a10 * B.a00 + A.a11 * B.a10 + A.a12 * B.a20;
    C.a11 = A.a10 * B.a01 + A.a11 * B.a11 + A.a12 * B.a21;
    C.a12 = A.a10 * B.a02 + A.a11 * B.a12 + A.a12 * B.a22;
    C.a20 = A.a20 * B.a00 + A.a21 * B.a10 + A.a22 * B.a20;
    C.a21 = A.a20 * B.a01 + A.a21 * B.a11 + A.a22 * B.a21;
    C.a22 = A.a20 * B.a02 + A.a21 * B.a12 + A.a22 * B.a22;
}
void MatrixTimes3313(const Matrix33 &A,const Matrix31 &B, Matrix31&C)
{
    C.x = A.a00 * B.x + A.a01 * B.y + A.a02 * B.z;
    C.y = A.a10 * B.x + A.a11 * B.y + A.a12 * B.z;
    C.z = A.a20 * B.x + A.a21 * B.y + A.a22 * B.z;
}

Matrix31 operator*(double k, const Matrix31 &v)
{
    return Matrix31(k * v.x, k * v.y, k * v.z);
}

Matrix31 operator*(const Matrix31 &v, double k)
{
    return Matrix31(k * v.x, k * v.y, k * v.z);
}
Matrix33 operator*(double k, const Matrix33 &M)
{
    return Matrix33(
        k * M.a00, k * M.a01, k * M.a02,
        k * M.a10, k * M.a11, k * M.a12,
        k * M.a20, k * M.a21, k * M.a22
    );
}

Matrix33 operator*(const Matrix33 &M, double k)
{
    return k * M; 
}
