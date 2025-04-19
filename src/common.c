#include "cstd.h"
float *ricker(int nt, float dt, float fm, float t0);
void pml_init();
void set_zero(float *p, int n);
void step_forward(float **p1, float **p2, float **vx1, float **vx2, float **vz1, float **vz2);
void write_float(const char *filename, float *p, int n);
static int nz, nx, nt;
static float dh, dt, fm, t0;
static float **pmlx, **pmlz;
int nb, nzb, nxb;
float v = 2500; /* velocity */
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
    if (!getparint("nb", &nb))
        nb = 30;
    t0 = 1. / fm;
    int sz = nz / 2 + nb, sx = nx / 2 + nb;
    nzb = nz + 2 * nb, nxb = nx + 2 * nb;
    /* wavelet */
    float *wt = ricker(nt, dt, fm, t0);
    float **p1 = alloc2float(nzb, nxb);
    float **p2 = alloc2float(nzb, nxb);
    float **vx1 = alloc2float(nzb, nxb);
    float **vx2 = alloc2float(nzb, nxb);
    float **vz1 = alloc2float(nzb, nxb);
    float **vz2 = alloc2float(nzb, nxb);
    float **tmp;
    set_zero(p1[0], nzb * nxb);
    set_zero(p2[0], nzb * nxb);
    set_zero(vx1[0], nzb * nxb);
    set_zero(vx2[0], nzb * nxb);
    set_zero(vz1[0], nzb * nxb);
    set_zero(vz2[0], nzb * nxb);
    pml_init();
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
    float px, pz;
    /* z-velocity */
    for (int ix = 0; ix < nzb; ix++)
    {
        for (int iz = 1; iz < nzb; iz++)
        {
            vz2[ix][iz] = ((p1[ix][iz] - p1[ix][iz - 1]) / dh - pmlz[ix][iz] * vz1[ix][iz]) * dt + vz1[ix][iz];
        }
    }
    /* x-velocity */
    for (int ix = 1; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nzb; iz++)
        {
            vx2[ix][iz] = ((p1[ix][iz] - p1[ix - 1][iz]) / dh - pmlx[ix][iz] * vx1[ix][iz]) * dt + vx1[ix][iz];
        }
    }
    /* stress */
    for (int ix = 1; ix < nxb - 1; ix++)
    {
        for (int iz = 1; iz < nzb - 1; iz++)
        {
            px = 0.5 * p1[ix][iz], pz = 0.5 * p1[ix][iz];
            px = (v * v * (vx2[ix + 1][iz] - vx2[ix][iz]) / dh - pmlx[ix][iz] * px) * dt + px;
            pz = (v * v * (vz2[ix][iz + 1] - vz2[ix][iz]) / dh - pmlz[ix][iz] * pz) * dt + pz;
            p2[ix][iz] = px + pz;
        }
    }
}
void write_float(const char *filename, float *p, int n)
{
    FILE *fp = efopen(filename, "wb");
    efwrite(p, sizeof(float), n, fp);
    efclose(fp);
}
void pml_init()
{
    const float D = nb;
    const float R = -1.5 * v / D * log(1e-5);
    pmlx = alloc2float(nzb, nxb);
    pmlz = alloc2float(nzb, nxb);
    for (int ix = 0; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nzb; iz++)
        {
            if (ix < nb)
                pmlx[ix][iz] = R * SQUARE((nb - ix) / D);
            else if (ix < nx + nb)
                pmlx[ix][iz] = 0;
            else
                pmlx[ix][iz] = R * SQUARE((ix + 1 - nb - nx) / D);
            if (iz < nb)
                pmlz[ix][iz] = R * SQUARE((nb - iz) / D);
            else if (iz < nz + nb)
                pmlz[ix][iz] = 0;
            else
                pmlz[ix][iz] = R * SQUARE((iz + 1 - nb - nz) / D);
        }
    }
}
void set_zero(float *p, int n)
{
    memset(p, 0, sizeof(float) * n);
}