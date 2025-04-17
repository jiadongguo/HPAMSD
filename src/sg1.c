/*
Fig. 1
*/
#include "cstd.h"
void get_sg_coeff2(int n, float *a);
void get_sg_coeff(int n, float *a);
float *ricker(int nt, float dt, float fm, float t0);
void step_forward(float **p0, float **p1, float **p2);
void write_float(const char *filename, float *p, int n);
float first_sg(float **p, int iz, int ix, int dir, int sign);
void laplace2(float **p, int ix, int iz, float *lap_x, float *lap_z);
static int order; /* differential accuracy L*/
static int nz, nx, nt;
static float dh, dt, fm, t0;
static float *c; /* first-order difference coefficient */
float v = 3000;  /* velocity */
int main(int argc, char **argv)
{
    initargs(argc, argv);
    char *wfd;
    if (!getparint("order", &order))
        order = 1;
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

    /* calculated difference coefficient */
    warn("Forward order=%d", order);
    order /= 2; /* order:2N */
    c = alloc1float(order);
    get_sg_coeff(order, c);
    /* wavelet */
    float *wt = ricker(nt, dt, fm, t0);
    float **p0 = alloc2float(nz, nx);
    float **p1 = alloc2float(nz, nx);
    float **p2 = alloc2float(nz, nx);
    float **tmp;
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
    float lap_x, lap_z;
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            laplace2(p1, ix, iz, &lap_x, &lap_z);
            p2[ix][iz] = SQUARE(v * dt / dh) * (lap_x + lap_z) + 2 * p1[ix][iz] - p0[ix][iz];
        }
    }
}
void laplace2(float **p, int ix, int iz, float *lap_x, float *lap_z)
{
    float pre, next;
    *lap_x = 0, *lap_z = 0;
    for (int i = 0; i < order; i++)
    {
        pre = first_sg(p, iz, ix + i, 2, 1);
        next = first_sg(p, iz, ix - i, 2, 0);
        *lap_x += c[i] * (pre - next);
        pre = first_sg(p, iz + i, ix, 1, 1);
        next = first_sg(p, iz - i, ix, 1, 0);
        *lap_z += c[i] * (pre - next);
    }
}
float first_sg(float **p, int iz, int ix, int dir, int sign)
{
    float ret = 0;
    int m;
    float pre, next;

    for (int i = 0; i < order; i++)
    {
        m = i + 1;
        pre = 0, next = 0;
        /* z */
        if (dir == 1)
        {
            if (sign == 1)
            {
                if (iz + m < nz && iz + m >= 0)
                    pre = p[ix][iz + m];
                if (iz - m + 1 >= 0 && iz - m + 1 < nz)
                    next = p[ix][iz - m + 1];
            }
            else
            {
                if (iz + m - 1 < nz && iz + m - 1 >= 0)
                    pre = p[ix][iz + m - 1];
                if (iz - m >= 0 && iz - m < nz)
                    next = p[ix][iz - m];
            }
        }
        /* x */
        else
        {

            if (sign == 1)
            {
                if (ix + m < nx && ix + m >= 0)
                    pre = p[ix + m][iz];
                if (ix - m + 1 >= 0 && ix - m + 1 < nx)
                    next = p[ix - m + 1][iz];
            }
            else
            {
                if (ix + m - 1 < nx && ix + m - 1 >= 0)
                    pre = p[ix + m - 1][iz];
                if (ix - m >= 0 && ix - m < nx)
                    next = p[ix - m][iz];
            }
        }
        ret += c[i] * (pre - next);
    }
    return ret;
}
void write_float(const char *filename, float *p, int n)
{
    FILE *fp = efopen(filename, "wb");
    efwrite(p, sizeof(float), n, fp);
    efclose(fp);
}