#ifndef __SIM_H__
#define __SIM_H__
typedef struct SIM
{
    int nz, nx;
    int nzb, nxb;
    int nbt, nbb, nbr, nbl;
    float dz, dx;
    int nt;
    float dt;
    float **v;
    float **vpad;
} SIM;
#endif