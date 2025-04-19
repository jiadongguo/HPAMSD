#include "cfft2.h"
#include "cstd.h"
int nz, nx;
float dh;
int sz, sx;
float c;
float dt;
float fm, t0;
int nt;
float *ricker(int nt, float dt, float fm, float t0);
void write_data(char *filename, float *d, int n);
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
    if (!getparfloat("c", &c))
        err("need c=");
    float *wt = ricker(nt, dt, fm, t0);
    float **p0 = alloc2float(nz, nx);
    float **p1 = alloc2float(nz, nx);
    float **p2 = alloc2float(nz, nx);
    float **tmp;
    memset(p0[0], 0, sizeof(float) * nz * nx);
    memset(p1[0], 0, sizeof(float) * nz * nx);
    memset(p2[0], 0, sizeof(float) * nz * nx);
    int nkz, nkx;
    cfft2_init(nz, nx, &nkz, &nkx);
    float *kx = alloc1float(nkx);
    float *kz = alloc1float(nkz);
    float kx0, kz0, dkx, dkz;
    kx0 = -0.5 / dh;
    kz0 = -0.5 / dh;
    dkx = 1 / (dh * nkx);
    dkz = 1 / (dh * nkz);
    for (int i = 0; i < nkx; i++)
        kx[i] = (kx0 + dkx * i) * 2 * PI;
    for (int i = 0; i < nkz; i++)
        kz[i] = (kz0 + dkz * i) * 2 * PI;
    fftwf_complex *kp1 = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * nkx * nkz);
    sz = nz / 2;
    sx = nx / 2;
    // p1[sx][sz] = 1;
    for (int it = 0; it < nt; it++)
    {
        if (it % 100 == 0)
        {
            printf("it=%d/%d\n", it, nt);
        }
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iz = 0; iz < nz; iz++)
            {
                kp1[ix * nz + iz][0] = p1[ix][iz];
                kp1[ix * nz + iz][1] = 0;
            }
        }
        cfft2(kp1, kp1);
        for (int ikx = 0; ikx < nkx; ikx++)
        {
            for (int ikz = 0; ikz < nkz; ikz++)
            {
                kp1[ikx * nkz + ikz][0] *= (kx[ikx] * kx[ikx] + kz[ikz] * kz[ikz]);
                kp1[ikx * nkz + ikz][1] *= (kx[ikx] * kx[ikx] + kz[ikz] * kz[ikz]);
            }
        }
        icfft2(kp1, kp1);
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iz = 0; iz < nz; iz++)
            {
                p2[ix][iz] = -(c * dt * c * dt) * kp1[ix * nz + iz][0] + 2 * p1[ix][iz] - p0[ix][iz];
            }
        }
        tmp = p0;
        p0 = p1;
        p1 = p2;
        p2 = tmp;
        p1[sx][sz] += wt[it];
    }
    write_data(wfd, p1[0], nz * nx);
    free(wt);
    free2(p0);
    free2(p1);
    free2(p2);
    return 0;
}
void write_data(char *filename, float *d, int n)
{
    FILE *fd = fopen(filename, "wb");
    fwrite(d, sizeof(float), n, fd);
    fclose(fd);
}
