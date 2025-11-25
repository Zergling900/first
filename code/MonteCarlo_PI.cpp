#include <random>
#include <iostream>
#include <cmath>

using namespace std;


int main()
{
    random_device z;
    mt19937 gen(z());
    uniform_real_distribution<> dis(0.0, 1.0);
    int N = 100000;
    int P = 0;

    for(int i = 0; i < N; i++)
    {
    double x = dis(gen);
    double y = dis(gen);
    if (x*x + y*y <= 1.0)
    {
        P++;
    }
    }

    double pi = 4.0 * P / N;
    cout << pi << endl;
    return 0;
}