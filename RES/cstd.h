#ifndef _cstd_h_
#define _cstd_h_
#include <unistd.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <cblas.h>
#ifndef NULL
#define NULL ((void *)0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef YES
#define YES (1)
#endif
#ifndef NO
#define NO (0)
#endif
#ifndef PI
#define PI 3.1415927
#endif

#ifndef PI2
#define PI2 9.8696044
#endif

#ifndef PI_1
#define PI_1 0.3183099
#endif

#ifndef PI_2
#define PI_2 0.1013212
#endif

#ifndef TOL
#define TOL 1e-12
#endif

#ifndef EPS
#define EPS FLT_EPSILON
#endif

#ifndef NULL
#define NULL ((void *)0)
#endif

#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef SQUARE
#define SQUARE(x) ((x) * (x))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef NINT
#define NINT(x) ((int)((x) > 0.0 ? (x) + 0.5 : (x) - 0.5))
#endif
#define STREQ(s, t) (strcmp(s, t) == 0)
#define STRLT(s, t) (strcmp(s, t) < 0)
#define STRGT(s, t) (strcmp(s, t) > 0)
#define DIM(a) (sizeof(a) / sizeof(a[0]))

#define iC(fun)            \
    {                      \
        int ierr = fun;    \
        assert(ierr == 0); \
    }
#define iA(expr)                                                             \
    {                                                                        \
        if ((expr) == 0)                                                     \
        {                                                                    \
            std::cerr << "wrong " << __LINE__ << " in " << __FILE__ << endl; \
            assert(expr);                                                    \
        }                                                                    \
    }
/* zpx data type */
typedef float cpx[2];
typedef double zpx[2];
/* interface of optimization algorithm */
typedef float (*opt_fg)(float *x, float *g);

/* IO interface */
void err(const char *fmt, ...);
void warn(const char *fmt, ...);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

/* allocate and free multi-dimensional arrays */
void *alloc1(size_t n1, size_t size);
void *realloc1(void *v, size_t n1, size_t size);
void **alloc2(size_t n1, size_t n2, size_t size);

void free1(void *p);
void free2(void **p);

int *alloc1int(size_t n1);
int *realloc1int(int *v, size_t n1);
int **alloc2int(size_t n1, size_t n2);
void free1int(int *p);
void free2int(int **p);

float *alloc1float(size_t n1);
float *realloc1float(float *v, size_t n1);
float **alloc2float(size_t n1, size_t n2);
void free1float(float *p);
void free2float(float **p);

double *alloc1double(size_t n1);
double *realloc1double(double *v, size_t n1);
double **alloc2double(size_t n1, size_t n2);
void free1double(double *p);
void free2double(double **p);

cpx *alloc1cpx(size_t n1);
cpx *realloc1cpx(cpx *v, size_t n1);
cpx **alloc2cpx(size_t n1, size_t n2);
void free1cpx(cpx *p);
void free2cpx(cpx **p);

zpx *alloc1zpx(size_t n1);
zpx *realloc1zpx(zpx *v, size_t n1);
zpx **alloc2zpx(size_t n1, size_t n2);
void free1zpx(zpx *p);
void free2zpx(zpx **p);

char *alloc1char(size_t n1);
char *realloc1char(char *v, size_t n1);
void free1char(char *p);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

/* system subroutine calls with error trapping */
FILE *efopen(const char *file, const char *mode);
FILE *efreopen(const char *file, const char *mode, FILE *stream1);
FILE *efdopen(int fd, const char *mode);
FILE *epopen(char *command, char *type);
int efclose(FILE *stream);
int epclose(FILE *stream);
int efflush(FILE *stream);
int eremove(const char *file);
int erename(const char *oldfile, const char *newfile);
int efseeko(FILE *stream, off_t offset, int origin);
int efseek(FILE *stream, off_t offset, int origin);
long eftell(FILE *stream);
off_t eftello(FILE *stream);
void erewind(FILE *stream);
int efgetpos(FILE *stream, fpos_t *position);
int efsetpos(FILE *stream, const fpos_t *position);
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream);
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream);

/* getpar parameter parsing */
void initargs(int argc, char **argv);
int getparint(char *name, int *p);
int getparuint(char *name, unsigned int *p);
int getparshort(char *name, short *p);
int getparushort(char *name, unsigned short *p);
int getparlong(char *name, long *p);
int getparulong(char *name, unsigned long *p);
int getparfloat(char *name, float *p);
int getpardouble(char *name, double *p);
int getparstring(char *name, char **p);
int getparstringarray(char *name, char **p);
int getnparint(int n, char *name, int *p);
int getnparuint(int n, char *name, unsigned int *p);
int getnparshort(int n, char *name, short *p);
int getnparushort(int n, char *name, unsigned short *p);
int getnparlong(int n, char *name, long *p);
int getnparulong(int n, char *name, unsigned long *p);
int getnparfloat(int n, char *name, float *p);
int getnpardouble(int n, char *name, double *p);
int getnparstring(int n, char *name, char **p);
int getnparstringarray(int n, char *name, char **p);
int getnpar(int n, char *name, char *type, void *ptr);
int countparname(char *name);
int countparval(char *name);
int countnparval(int n, char *name);
void checkpars(void);

#endif
