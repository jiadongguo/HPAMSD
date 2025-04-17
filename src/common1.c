/*
普通二阶中心差分
*/
#include "cstd.h"
float *ricker(int nt, float dt, float fm, float t0);
void set_zero(float *p, int n);
void step_forward(float **p0, float **p1, float **p2);
void write_float(const char *filename, float *p, int n);
float bnd(float **p, int ix, int iz, int nx, int nz);
static int nz, nx, nt;
static float dh, dt, fm, t0;
float v = 3000; /* velocity */
int main(int argc, char **argv)
{
    initargs(argc, argv);
    char *wfd;
    if (!getparint("nz", &nz))
        err("need nz=");
    if (!getparint("nx", &nx))
        err("need nx=");
    if (!getparint("nt", &nt))
        err("need nt=");
    if (!getparfloat("dh", &dh))
        err("need dh=");
    if (!getparfloat("dt", &dt))
        err("dt=");
    if (!getparfloat("fm", &fm))
        err("need fm=");
    if (!getparfloat("t0", &t0))
        t0 = 1. / fm;
    if (!getparstring("wfd", &wfd))
        err("need wfd=");
    int sz = nz / 2, sx = nx / 2;
    /* wavelet */
    float *wt = ricker(nt, dt, fm, t0);
    float **p0 = alloc2float(nz, nx);
    float **p1 = alloc2float(nz, nx);
    float **p2 = alloc2float(nz, nx);
    float **tmp;
    set_zero(p0[0], nz * nx);
    set_zero(p1[0], nz * nx);
    set_zero(p2[0], nz * nx);
    for (int it = 0; it < nt; it++)
    {
        if (it % 100 == 0)
        {
            warn("it=%d/%d", it, nt);
        }
        step_forward(p0, p1, p2);
        tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
        p1[sx][sz] += wt[it];
    }
    write_float(wfd, p1[0], nz * nx);
    return 0;
}
void step_forward(float **p0, float **p1, float **p2)
{
    /* stress */
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            p2[ix][iz] = SQUARE(v * dt / dh) *
                             ((bnd(p1, ix + 1, iz, nx, nz) - 2 * p1[ix][iz] + bnd(p1, ix - 1, iz, nx, nz)) +
                              (bnd(p1, ix, iz + 1, nx, nz) - 2 * p1[ix][iz] + bnd(p1, ix, iz - 1, nx, nz))) +
                         2 * p1[ix][iz] - p0[ix][iz];
        }
    }
}
void write_float(const char *filename, float *p, int n)
{
    FILE *fp = efopen(filename, "wb");
    efwrite(p, sizeof(float), n, fp);
    efclose(fp);
}
void set_zero(float *p, int n)
{
    memset(p, 0, sizeof(float) * n);
}
float bnd(float **p, int ix, int iz, int nx, int nz)
{
    if (ix < 0 || ix >= nx || iz < 0 || iz >= nz)
        return 0;
    return p[ix][iz];
}