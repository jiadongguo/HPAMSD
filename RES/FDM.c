/*
普通中心差分
*/
#include "cstd.h"
void set_zero(float *p, int n);
void step_forward(float **p0, float **p1, float **p2);
void write_float(const char *filename, float *p, int n);
float bnd(float **p, int ix, int iz);
float *ricker(int nt, float dt, float fm, float t0);
float laplace2(float **p, int ix, int iz);
static int nz, nx, nt;
static float dh, dt, fm, t0;
static int order;
float v; /* velocity */
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
    if (!getparfloat("c", &v))
        err("need c=");
    if (!getparint("order", &order))
        order = 10;
    warn("order=%d", order);
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
    free2(p0);
    free2(p1);
    free2(p2);
    free(wt);
    return 0;
}
void step_forward(float **p0, float **p1, float **p2)
{
    /* stress */
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            p2[ix][iz] = SQUARE(v * dt / dh) * laplace2(p1, ix, iz) + 2 * p1[ix][iz] - p0[ix][iz];
            // p2[ix][iz] = SQUARE(v * dt / dh) * (p1[ix + 1][iz] + p1[ix - 1][iz] + p1[ix][iz + 1] + p1[ix][iz - 1] - 4 * p1[ix][iz]) + 2 * p1[ix][iz] - p0[ix][iz];
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
float bnd(float **p, int ix, int iz)
{
    if (ix < 0 || ix >= nx || iz < 0 || iz >= nz)
        return 0;
    return p[ix][iz];
}
float laplace2(float **p, int ix, int iz)
{
    float *c;
    static float c2[] = {-2., 1.};
    static float c4[] = {-5. / 2, 4. / 3, -1. / 12};
    static float c6[] = {-49. / 18, 3. / 2, -3. / 20, 1. / 90};
    static float c8[] = {-205. / 72, 8. / 5, -1. / 5, 8. / 315, -1. / 560};
    static float c10[] = {-5269. / 1800, 5. / 3, -5. / 21, 5. / 126, -5. / 1008, 1. / 3150};
    switch (order)
    {
    case 2:
        c = c2;
        break;
    case 4:
        c = c4;
        break;
    case 6:
        c = c6;
        break;
    case 8:
        c = c8;
    default:
        c = c10;
        break;
    }
    float ret = 2 * c[0] * p[ix][iz];
    for (int i = 1; i < order / 2 + 1; i++)
    {
        ret += c[i] * (bnd(p, ix + i, iz) + bnd(p, ix - i, iz) + bnd(p, ix, iz + i) + bnd(p, ix, iz - i));
    }
    return ret;
}