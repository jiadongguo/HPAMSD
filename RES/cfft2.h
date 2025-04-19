#ifndef __CFFT2_H__
#define __CFFT2_H__
#include <fftw3.h>
#include "cstd.h"
void cfft2_init(int n1, int n2, int *nk1, int *nk2);
void cfft2(fftwf_complex *st, fftwf_complex *sf);
void icfft2(fftwf_complex *sf, fftwf_complex *st);
void cfft_close();
#endif