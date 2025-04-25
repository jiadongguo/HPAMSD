#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define EPS 1E-7
#define pi 3.1415926
#define SQUARE(x) ((x) * (x))
void mkmodel();
float *allocfloat(int n);
float **alloc2float(int n1, int n2);
float *ricker(int nt, float dt, float fm, float t0);
void seis_record(float **p, float *rcd);
void pml_init();
void step_forward(float **p0, float **p1, float **p2);
void add_source(int is, float **p, float *src, int flag);
void corr(float **img, float **p1, float **p2);
/*< grid size >*/
int nb = 30;
int nx = 400, nz = 400;
float dx = 10, dz = 10;
int nxb, nzb;
/*< source position >*/
int sz = 10, sx = 10, jsx = 20;
int ns = 15;
/*< receiver position >*/
int rz = 5, rx = 0, jrx = 1;
int nr = 400;
/*< simulation time >*/
int nt = 500;
float dt = 0.001;
/*< wavelet  parameters>*/
float fm = 20;
float t0 = 0.1; // delay
/*< velocity model >*/
float **c;
float **rcd;
/*< PML parameter>*/
float **d;
/*< imag >*/
float **img;
int main()
{
    nzb = nz + 2 * nb;
    nxb = nx + 2 * nb;
    img = alloc2float(nzb, nxb);
    float **tmpimg = alloc2float(nzb, nxb);
    float **illum = alloc2float(nzb, nxb);
    mkmodel();
    pml_init();
    rcd = alloc2float(nr, nt);
    float *wt = ricker(nt, dt, fm, t0);
    float **p0 = alloc2float(nzb, nxb);
    float **p1 = alloc2float(nzb, nxb);
    float **p2 = alloc2float(nzb, nxb);
    float **tmp;
    FILE *fp;
    FILE *frcd = fopen("rcd", "wb");
    memset(img[0], 0, sizeof(float) * nzb * nxb);
    for (int is = 0; is < ns; is++)
    {
        memset(p0[0], 0, sizeof(float) * nzb * nxb);
        memset(p1[0], 0, sizeof(float) * nzb * nxb);
        memset(p2[0], 0, sizeof(float) * nzb * nxb);
        memset(illum[0], 0, sizeof(float) * nzb * nxb);
        /* forward */
        fp = fopen("seis", "wb");
        for (int it = 0; it < nt; it++)
        {
            if (it % 100 == 0)
            {
                printf("Forward is=%d/%d,it=%d/%d\n", is, ns, it, nt);
            }
            step_forward(p0, p1, p2);
            tmp = p0;
            p0 = p1;
            p1 = p2;
            p2 = tmp;
            add_source(is, p1, &wt[it], 1);
            seis_record(p1, rcd[it]);
            fwrite(p1[0], sizeof(float), nzb * nxb, fp);
            for (int ix = 0; ix < nzb; ix++)
            {
                for (int iz = 0; iz < nzb; iz++)
                {
                    illum[ix][iz] += SQUARE(p1[ix][iz]);
                }
            }
        }
        fwrite(rcd[0], sizeof(float), nr * nt, frcd);
        fclose(fp);
        /* backward */
        memset(p0[0], 0, sizeof(float) * nzb * nxb);
        memset(p1[0], 0, sizeof(float) * nzb * nxb);
        memset(p2[0], 0, sizeof(float) * nzb * nxb);
        fp = fopen("seis", "rb");
        for (int it = nt - 1; it >= 0; it--)
        {
            if (it % 100 == 0)
            {
                printf("Backward is=%d/%d,it=%d/%d\n", is, ns, it, nt);
            }
            fseek(fp, it * sizeof(float) * nzb * nxb, SEEK_SET);
            step_forward(p0, p1, p2);
            tmp = p0;
            p0 = p1;
            p1 = p2;
            p2 = tmp;
            add_source(is, p1, rcd[it], 0);
            fread(p2[0], sizeof(float), nzb * nxb, fp);
            corr(tmpimg, p2, p1);
            for (int ix = 0; ix < nxb; ix++)
            {
                for (int iz = 0; iz < nzb; iz++)
                {
                    img[ix][iz] += tmpimg[ix][iz] ;// (illum[ix][iz] + EPS);
                }
            }
        }

        fclose(fp);
    }
    fp = fopen("img", "wb");
    fwrite(img[0], sizeof(float), nzb * nxb, fp);
    fclose(fp);
    fclose(frcd);
    return 0;
}
float *allocfloat(int n)
{
    float *p = (float *)malloc(sizeof(float) * n);
    return p;
}
float **alloc2float(int n1, int n2)
{
    float **p = (float **)malloc(sizeof(float *) * n2);
    p[0] = (float *)malloc(sizeof(float) * n1 * n2);
    for (int i = 1; i < n2; i++)
    {
        p[i] = p[0] + i * n1;
    }

    return p;
}
void mkmodel()
{
    c = alloc2float(nzb, nxb);
    for (int ix = 0; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nzb; iz++)
        {
            if (iz < nzb / 2)
                c[ix][iz] = 2500;
            else
                c[ix][iz] = 3000;
        }
    }
}
float *ricker(int nt, float dt, float fm, float t0)
{
    float *wt = allocfloat(nt);
    float tmp;
    for (int it = 0; it < nt; it++)
    {
        tmp = fm * pi * (it * dt - t0);
        tmp *= tmp;
        wt[it] = (1 - 2 * tmp) * exp(-tmp);
    }
    return wt;
}
void pml_init()
{
    d = alloc2float(nzb, nxb);
    // memset(d[0], 0, sizeof(float) * nzb * nxb);
    // return;
    float ref = 1e-6;
    float alpha = 1.5 * log(1 / ref);
    float width, dis;
    /* inner */
    for (int ix = nb; ix < nb + nx; ix++)
    {
        for (int iz = nb; iz < nb + nz; iz++)
        {
            d[ix][iz] = 0;
        }
    }

    width = nb * dz;
    for (int ix = nb; ix < nb + nx; ix++)
    {
        /*top */
        for (int iz = 0; iz < nb; iz++)
        {
            dis = (nb - iz) * dz;
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
        /*bottom*/
        for (int iz = nb + nz; iz < nzb; iz++)
        {
            dis = (iz + 1 - nb - nz) * dz;
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
    }
    width = nb * dx;
    for (int iz = nb; iz < nb + nz; iz++)
    {
        /* left */
        for (int ix = 0; ix < nb; ix++)
        {
            dis = (nb - ix) * dx;
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
        /* right */
        for (int ix = nb + nx; ix < nxb; ix++)
        {
            dis = (ix + 1 - nb - nx) * dx;
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
    }
    /* corner */
    width = hypot(nb * dx, nb * dz);
    for (int ix = 0; ix < nb; ix++)
    {
        for (int iz = 0; iz < nb; iz++)
        {
            dis = hypot((nb - ix) * dx, (nb - iz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
        for (int iz = nb + nz; iz < nzb; iz++)
        {
            dis = hypot((nb - ix) * dx, (iz + 1 - nb - nz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
    }
    for (int ix = nb + nx; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nb; iz++)
        {
            dis = hypot((ix + 1 - nb - nx) * dx, (nb - iz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
        for (int iz = nb + nz; iz < nzb; iz++)
        {
            dis = hypot((ix + 1 - nb - nx) * dx, (iz + 1 - nb - nz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * c[ix][iz];
        }
    }
}
void step_forward(float **p0, float **p1, float **p2)
{
    float dxx, dzz;
    float A, M, N;
    for (int ix = 1; ix < nxb - 1; ix++)
    {
        for (int iz = 1; iz < nzb - 1; iz++)
        {
            dxx = (p1[ix + 1][iz] - 2 * p1[ix][iz] + p1[ix - 1][iz]) / (dx * dx);
            dzz = (p1[ix][iz + 1] - 2 * p1[ix][iz] + p1[ix][iz - 1]) / (dz * dz);
            A = 1 + d[ix][iz] * dt;
            M = 2 - d[ix][iz] * d[ix][iz] * dt * dt;
            N = d[ix][iz] * dt - 1;
            p2[ix][iz] = M / A * p1[ix][iz] + N / A * p0[ix][iz] + c[ix][iz] * c[ix][iz] * dt * dt * (dxx + dzz);
        }
    }
}
void seis_record(float **p, float *rcd)
{
    for (int ir = 0; ir < nr; ir++)
    {
        rcd[ir] = p[rx + nb + ir * jrx][rz + nb];
    }
}
/* <flag=1:forward,flag=0:adjoint >*/
void add_source(int is, float **p, float *src, int flag)
{
    if (flag)
    {
        p[nb + sx + is * jsx][nb + sz] += src[0];
    }
    else
    {
        for (int ir = 0; ir < nr; ir++)
        {
            p[rx + nb + ir * jrx][rz + nb] += src[ir];
        }
    }
}
void corr(float **img, float **p1, float **p2)
{
    for (int ix = 0; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nzb; iz++)
        {
            img[ix][iz] = p1[ix][iz] * p2[ix][iz];
        }
    }
}