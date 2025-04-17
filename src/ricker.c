#include "cstd.h"
float *ricker(int nt, float dt, float fm, float t0)
{
    float *wt = alloc1float(nt);
    float tmp;
    for (int it = 0; it < nt; it++)
    {
        tmp = fm * PI * (it * dt - t0);
        tmp *= tmp;
        wt[it] = (1 - 2 * tmp) * exp(-tmp);
    }
    return wt;
}