#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
static int pw2[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
static int nx, nz;
static int nxz;
static int nkx, nkz;
static int nkxz;
static fftwf_complex *d;
fftwf_plan forward, backward;
static int fft_next_fast_size(int n)
{
    int j = 0;
    while (pw2[j] < n)
    {
        j++;
    }
    return pw2[j];
}
void cfft2_init(int n1, int n2, int *nk1, int *nk2)
{
    nz = n1;
    nx = n2;
    nkx = fft_next_fast_size(n2);
    nkz = fft_next_fast_size(n1);
    nkxz = nkx * nkz;
    d = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * nkxz);
    *nk1 = nkz;
    *nk2 = nkx;
    forward = fftwf_plan_dft_2d(nkx, nkz, d, d, FFTW_FORWARD, FFTW_MEASURE);
    backward = fftwf_plan_dft_2d(nkx, nkz, d, d, FFTW_BACKWARD, FFTW_MEASURE);
}
void cfft2(fftwf_complex *st /*[nx][nz] */, fftwf_complex *sf /*[nkx][nkz]*/)
{
    memset(d, 0, sizeof(fftwf_complex) * nkxz);
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            if ((i + j) % 2 == 0)
            {
                d[i * nkz + j][0] = st[i * nz + j][0];
                d[i * nkz + j][1] = st[i * nz + j][1];
            }
            else
            {
                d[i * nkz + j][0] = -st[i * nz + j][0];
                d[i * nkz + j][1] = -st[i * nz + j][1];
            }
        }
    }
    fftwf_execute(forward);
    memcpy(sf, d, sizeof(fftwf_complex) * nkxz);
}
void icfft2(fftwf_complex *sf /*[nkx][nkz]*/, fftwf_complex *st /*[nx][nz] */)
{
    memcpy(d, sf, sizeof(fftwf_complex) * nkxz);
    fftwf_execute(backward);
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            if ((i + j) % 2 == 0)
            {
                st[i * nz + j][0] = d[i * nkz + j][0] / nkxz;
                st[i * nz + j][1] = d[i * nkz + j][1] / nkxz;
            }
            else
            {
                st[i * nz + j][0] = -d[i * nkz + j][0] / nkxz;
                st[i * nz + j][1] = -d[i * nkz + j][1] / nkxz;
            }
        }
    }
}
void cfft_close()
{
    fftwf_free(d);
    fftwf_destroy_plan(forward);
    fftwf_destroy_plan(backward);
}