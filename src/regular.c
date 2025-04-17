/*
Figure 3:central difference PML
*/
#include "cstd.h"
void get_reg_coeff(int n, float *a);
float *ricker(int nt, float dt, float fm, float t0);
void window2d(float **a, float **b);
void expand2d(float **b, float **a);
float first_order(float **p, int ix, int iz, int dir);
void step_forward(float **p1, float **p2, float **vx1, float **vx2, float **vz1, float **vz2);
void write_float(const char *filename, float *p, int n);
static int order; /* differential accuracy L*/
static int nz, nx, nt;
static int nzb, nxb;
static int nb;
static float dh, dt, fm, t0;
static float *c; /* first-order difference coefficient */
float v = 2500;  /* velocity */
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
    if (!getparint("nb", &nb))
        nb = 30;
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
    int sz = nz / 2 + nb, sx = nx / 2 + nb;

    /* calculated difference coefficient */
    order /= 2;
    c = alloc1float(order);
    get_reg_coeff(order, c);
    /* wavelet */
    float *wt = ricker(nt, dt, fm, t0);
    /* pad grid size */
    nzb = nz + 2 * nb;
    nxb = nx + 2 * nb;
    float **p1 = alloc2float(nzb, nxb);
    float **p2 = alloc2float(nzb, nxb);
    float **vx1 = alloc2float(nzb, nxb);
    float **vx2 = alloc2float(nzb, nxb);
    float **vz1 = alloc2float(nzb, nxb);
    float **vz2 = alloc2float(nzb, nxb);
    float **tmp;
    for (int it = 0; it < nt; it++)
    {
        if (it % 100 == 0)
        {
            warn("it=%d/%d", it, nt);
        }
        step_forward(p1, p2, vx1, vx2, vz1, vz2);
        tmp = p1, p1 = p2, p2 = tmp;
        tmp = vx1, vx1 = vx2, vx2 = tmp;
        tmp = vz1, vz1 = vz2, vz2 = tmp;
        p1[sx][sz] += wt[it];
    }
    write_float(wfd, p1[0], nzb * nxb);
    return 0;
}
void step_forward(float **p1, float **p2, float **vx1, float **vx2, float **vz1, float **vz2)
{
    float delta_x, delta_z;
    const float R = (3 * v) / (2 * nb * dh) * log(1 / (1e-7));
    const float w = 1. * nb;
    float dpx, dpz;
    float px, pz;
    for (int ix = 0; ix < nxb; ix++)
    {
        /* determine the position in the x direction */
        if (ix < nb)
            delta_x = R * SQUARE((nb - ix) / w);
        else if (ix < nb + nx)
            delta_x = 0;
        else
            delta_x = R * SQUARE((ix - nb - nx + 1) / w);
        for (int iz = 0; iz < nzb; iz++)
        {
            /* determine the position in the z direction */
            if (iz < nb)
                delta_z = R * SQUARE((nb - iz) / w);
            else if (iz < nb + nz)
                delta_z = 0;
            else
                delta_z = R * SQUARE((iz - nb - nz + 1) / w);
            dpx = first_order(p1, ix, iz, 2);
            dpz = first_order(p1, ix, iz, 1);
            /* split */
            px = 0.5 * p1[ix][iz];
            pz = 0.5 * p1[ix][iz];
            /* velocity px=pz=0.5*p */
            vx2[ix][iz] = (dpx / dh - delta_x * vx1[ix][iz]) * dt + vx1[ix][iz];
            vz2[ix][iz] = (dpz / dh - delta_z * vz1[ix][iz]) * dt + vz1[ix][iz];
            /* stress */
            px = (first_order(vx1, ix, iz, 2) / dh * v * v - delta_x * px) * dt + px;
            pz = (first_order(vz1, ix, iz, 1) / dh * v * v - delta_z * pz) * dt + pz;
            p2[ix][iz] = px + pz;
        }
    }
}
void window2d(float **a, float **b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
    int iz, ix;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix, iz) \
    shared(b, a, nb, nz, nx)
#endif
    for (ix = 0; ix < nx; ix++)
    {
        for (iz = 0; iz < nz; iz++)
        {
            a[ix][iz] = b[nb + ix][nb + iz];
        }
    }
}
void expand2d(float **b, float **a)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int iz, ix;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix, iz) \
    shared(b, a, nb, nz, nx)
#endif
    for (ix = 0; ix < nx; ix++)
    {
        for (iz = 0; iz < nz; iz++)
        {
            b[nb + ix][nb + iz] = a[ix][iz];
        }
    }

    for (ix = 0; ix < nxb; ix++)
    {
        for (iz = 0; iz < nb; iz++)
        {
            b[ix][iz] = b[ix][nb];
            b[ix][nzb - iz - 1] = b[ix][nzb - nb - 1];
        }
    }

    for (ix = 0; ix < nb; ix++)
    {
        for (iz = 0; iz < nzb; iz++)
        {
            b[ix][iz] = b[nb][iz];
            b[nxb - ix - 1][iz] = b[nxb - nb - 1][iz];
        }
    }
}

/* calculate the first-order difference, it is necessary to divide by the interval */
float first_order(float **p, int ix, int iz, int dir /* 1:z-direction 2:x-direction */)
{
    int m;
    float ret = 0, k1, k2;
    // if (dir < 1 || dir > 2)
    //     err("dir must be 1 or 2 in function first_order");
    for (int i = 0; i < order; i++)
    {
        m = i + 1;
        if (dir == 1)
        {
            if (iz - m < 0)
                k1 = 0;
            else
                k1 = p[ix][iz - m];
            if (iz + m >= nzb)

                k2 = 0;
            else
                k2 = p[ix][iz + m];
        }
        else
        {
            if (ix - m < 0)
                k1 = 0;
            else
                k1 = p[ix - m][iz];
            if (ix + m >= nxb)

                k2 = 0;
            else
                k2 = p[ix + m][iz];
        }
        ret += c[i] * (k2 - k1);
    }
    return ret;
}
void write_float(const char *filename, float *p, int n)
{
    FILE *fp = efopen(filename, "wb");
    efwrite(p, sizeof(float), n, fp);
    efclose(fp);
}