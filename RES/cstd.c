#include "cstd.h"
static char *prog = NULL;
void err(const char *fmt, ...)
/*< Outputs an error message to stderr and terminates the program.
---
Format and variable arguments follow printf convention. Additionally, a ':' at
the end of format adds system information for system errors. >*/
{
    va_list args;
    int last;
    (void)fflush(stdout);
    va_start(args, fmt);

    last = strlen(fmt) - 1;

    /* print out name of program causing error */
    fprintf(stderr, "%s: ", prog);

    /* print out remainder of message */
    (void)vfprintf(stderr, fmt, args);
    va_end(args);

    fprintf(stderr, "\n");
    (void)fflush(stderr);

    exit(EXIT_FAILURE);
}
void warn(const char *fmt, ...)
/* Outputs a warning message to stderr.
---
Format and variable arguments follow printf convention. Additionally, a ':' at
the end of format adds system information for system errors. */
{
    va_list args;
    int last;
    (void)fflush(stdout);
    va_start(args, fmt);
    last = strlen(fmt) - 1;

    /* if fmt ends with ';', carriage return */
    if (fmt[0] != '\0' && fmt[last] == ';')
        fprintf(stderr, "\r");

    /* if fmt starts with '.', new line */
    if (fmt[0] == '.')
    {
        fprintf(stderr, "\n");
        if (0 == last)
        {
            (void)fflush(stderr);
            return;
        }
    }

    /* print out name of program causing error */
    fprintf(stderr, "%s: ", prog);

    /* print out remainder of message */
    (void)vfprintf(stderr, fmt, args);
    va_end(args);

    /* if fmt ends with ';', do not end line */
    if (fmt[0] == '\0' || fmt[last] != ';')
        fprintf(stderr, "\n");

    (void)fflush(stderr);
}
/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************
ATOPKGE - convert ascii to arithmetic and with error checking


eatoh		ascii to short
eatou		ascii to unsigned short
eatoi		ascii to int
eatop		ascii to unsigned
eatol		ascii to long
eatov		ascii to unsigned long
eatof		ascii to float
eatod		ascii to double

****************************************************************************
Function Prototypes:
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

****************************************************************************
Input:
s		string

Returned:	type indicated

****************************************************************************
Notes:
Each of these routines acts like atoi, but has error checking:

This is a major revision of the tedious code required before
vendors implemented the ANSI C strtol, strtoul and strtod.

In addition to the size checks for each integer type, a
specific test on errno is required.  For example, INT_MAX
may (and probably does) equal LONG_MAX.  In this case,
if fed a number exceeding INT_MAX (and LONG_MAX), strtol
will make a quiet return with the wrong answer and it is up
to the user to check if errno == ERANGE.

Size limits are machine dependent and are read from the
ANSI C include files limits.h and float.h.

Bug Report: With NeXT c and Gnucc, when x > DBL_MAX (x <-DBL_MAX),
the return value from strtod was +Infinity (-Infinity), not HUGE_VAL
and more important, errno was not set to ERANGE.  To cope with this,
I put explicit size checks in eatod (which would not be needed if
errno were set as it should be in ANSI C.    jkc 01/29/94

On IBM RS6000, the return value from strtod was +-Inf on
overflow, but errno was set correctly.

****************************************************************************
References:
For old code:
Plum: Reliable Data Structures in C, p. 2-17.
Kernighan and Ritchie: The C Programming Language, p. 58.

CWP: Jack K. Cohen, Brian Sumner

For new code:
ANSI C routines with a little help from Jack

****************************************************************************
Author: Jack Cohen, Center for Wave Phenomena, 1994.
***************************************************************************/
/**************** end self doc ********************************/

/* eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX} */
short eatoh(char *s)
{
    long n = strtol(s, NULL, 10);

    if ((n > SHRT_MAX) || (n < SHRT_MIN) || (errno == ERANGE))
        err("%s: eatoh: overflow", __FILE__);

    return (short)n;
}
/* eatou - convert string s to unsigned short integer {0:USHRT_MAX} */
unsigned short eatou(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);

    if ((n > USHRT_MAX) || (errno == ERANGE))
        err("%s: eatou: overflow", __FILE__);

    return (unsigned short)n;
} /* eatoi - convert string s to integer {INT_MIN:INT_MAX} */
int eatoi(char *s)
{
    long n = strtol(s, NULL, 10);

    if ((n > INT_MAX) || (n < INT_MIN) || (errno == ERANGE))
        err("%s: eatoi: overflow", __FILE__);

    return (int)n;
}

/* eatop - convert string s to unsigned integer {0:UINT_MAX} */
unsigned int eatop(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);

    if ((n > UINT_MAX) || (errno == ERANGE))
        err("%s: eatop: overflow", __FILE__);

    return (unsigned int)n;
}

/* eatol - convert string s to long integer {LONG_MIN:LONG_MAX} */
long eatol(char *s)
{
    long n = strtol(s, NULL, 10);

    if (errno == ERANGE)
        err("%s: eatol: overflow", __FILE__);

    return n;
}

/* eatov - convert string s to unsigned long {0:ULONG_MAX} */
unsigned long eatov(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);

    if (errno == ERANGE)
        err("%s: eatov: overflow", __FILE__);

    return n;
}

/* eatof - convert string s to float {-FLT_MAX:FLT_MAX} */
float eatof(char *s)
{
    float x = strtod(s, NULL);

    if ((x > FLT_MAX) || (x < -FLT_MAX) || (errno == ERANGE))
        err("%s: eatof: overflow", __FILE__);

    return (float)x;
}

/* eatod - convert string s to double {-DBL_MAX:DBL_MAX} */
double eatod(char *s)
{
    double x = strtod(s, NULL);

    /* errno == ERANGE suffices if compiler sets errno on overflow */
    if ((errno == ERANGE) || (x > DBL_MAX) || (x < -DBL_MAX))
        err("%s: eatod: overflow", __FILE__);

    return x;
}
/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
ALLOC - Allocate and free multi-dimensional arrays

alloc1		allocate a 1-d array
realloc1	re-allocate a 1-d array
free1		free a 1-d array
alloc2		allocate a 2-d array
free2		free a 2-d array
alloc3		allocate a 3-d array
free3		free a 3-d array
alloc4		allocate a 4-d array
free4		free a 4-d array
alloc5		allocate a 5-d array
free5		free a 5-d array
alloc6		allocate a 6-d array
free6		free a 6-d arrayalloc1int
allocate a 1-d array of ints
realloc1int	re-allocate a 1-d array of ints
free1int	free a 1-d array of ints
alloc2int	allocate a 2-d array of ints
free2int	free a 2-d array of ints
alloc3int	allocate a 3-d array of ints
free3int	free a 3-d array of ints
alloc1float	allocate a 1-d array of floats
realloc1float	re-allocate a 1-d array of floats
free1float	free a 1-d array of floats
alloc2float	allocate a 2-d array of floats
free2float	free a 2-d array of floats
alloc3float	allocate a 3-d array of floats
free3float	free a 3-d array of floats
alloc4float	allocate a 4-d array of floats
free4float      free a 4-d array of floats
alloc5float     allocate a 5-d array of floats
free5float      free a 5-d array of floats
alloc6float     allocate a 6-d array of floats
free6float      free a 6-d array of floats
alloc4int       allocate a 4-d array of ints
free4int        free a 4-d array of ints
alloc5int       allocate a 5-d array of ints
free5int        free a 5-d array of ints
alloc5uchar	allocate a 5-d array of unsigned chars
free5uchar	free a 5-d array of unsiged chars
alloc5ushort    allocate a 5-d array of unsigned shorts
free5ushort     free a 5-d array of unsiged shorts
alloc6ushort    allocate a 6-d array of unsigned shorts
free6ushort     free a 6-d array of unsiged shorts
alloc1double	allocate a 1-d array of doubles
realloc1double	re-allocate a 1-d array of doubles
free1double	free a 1-d array of doubles
alloc2double	allocate a 2-d array of doubles
free2double	free a 2-d array of doubles
alloc3double	allocate a 3-d array of doubles
free3double	free a 3-d array of doubles
alloc1cpx	allocate a 1-d array of cpxs
realloc1cpx	re-allocate a 1-d array of cpxs
free1cpx	free a 1-d array of cpxs
alloc2cpx	allocate a 2-d array of cpxs
free2cpx	free a 2-d array of cpxs
alloc3cpx	allocate a 3-d array of cpxs
free3cpx	free a 3-d array of cpxs

alloc1zpx   allocate a 1-d array of cpxs
realloc1zpx re-allocate a 1-d array of cpxs
free1zpx    free a 1-d array of cpxs
alloc2zpx   allocate a 2-d array of cpxs
free2zpx    free a 2-d array of cpxs
alloc3zpx   allocate a 3-d array of cpxs
free3zpx    free a 3-d array of cpxs

******************************************************************************
Function Prototypes:
void *alloc1 (size_t n1, size_t size);
void *realloc1 (void *v, size_t n1, size_t size);
void free1 (void *p);
void **alloc2 (size_t n1, size_t n2, size_t size);
void free2 (void **p);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void free3 (void ***p);
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void free4 (void ****p);
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void free5 (void *****p);
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                  size_t size);
void free6 (void ******p);
int *alloc1int (size_t n1);
int *realloc1int (int *v, size_t n1);
void free1int (int *p);
int **alloc2int (size_t n1, size_t n2);
void free2int (int **p);
int ***alloc3int (size_t n1, size_t n2, size_t n3);
void free3int (int ***p);
float *alloc1float (size_t n1);
float *realloc1float (float *v, size_t n1);
void free1float (float *p);
float **alloc2float (size_t n1, size_t n2);
void free2float (float **p);
float ***alloc3float (size_t n1, size_t n2, size_t n3);
void free3float (float ***p);
float ****alloc4float (size_t n1, size_t n2, size_t n3, size_t n4);
void free4float (float ****p);
float *****alloc5float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5float (float *****p);
float ******alloc6float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5,
                         size_t n6);
void free6float (float ******p);
int ****alloc4int (size_t n1, size_t n2, size_t n3, size_t n4);
void free4int (int ****p);
int *****alloc5int (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5int (int *****p);
unsigned char *****alloc5uchar(size_t n1,size_t n2,size_t n3,size_t n4,
    size_t n5);
void free5uchar(unsigned char *****p);
unsigned short *****alloc5ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5);
void free5ushort(unsigned short *****p);
unsigned short ******alloc6ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5,size_t n6);
void free6ushort(unsigned short ******p);
double *alloc1double (size_t n1);
double *realloc1double (double *v, size_t n1);
void free1double (double *p);
double **alloc2double (size_t n1, size_t n2);
void free2double (double **p);
double ***alloc3double (size_t n1, size_t n2, size_t n3);
void free3double (double ***p);
cpx *alloc1cpx (size_t n1);
cpx *realloc1cpx (cpx *v, size_t n1);
void free1cpx (cpx *p);
cpx **alloc2cpx (size_t n1, size_t n2);
void free2cpx (cpx **p);
cpx ***alloc3cpx (size_t n1, size_t n2, size_t n3);
void free3cpx (cpx ***p);

cpx *alloc1zpx (size_t n1);
cpx *realloc1zpx (zpx *v, size_t n1);
void free1zpx (zpx *p);
cpx **alloc2zpx (size_t n1, size_t n2);
void free2zpx (zpx **p);
cpx ***alloc3zpx (size_t n1, size_t n2, size_t n3);
void free3zpx (zpx ***p);


******************************************************************************
Notes:
The functions defined below are intended to simplify manipulation
of multi-dimensional arrays in scientific programming in C.  These
functions are useful only because true multi-dimensional arrays
in C cannot have variable dimensions (as in FORTRAN).  For example,
the following function IS NOT valid in C:
    void badFunc(a,n1,n2)
    float a[n2][n1];
    {
        a[n2-1][n1-1] = 1.0;
    }
However, the following function IS valid in C:
    void goodFunc(a,n1,n2)
    float **a;
    {
        a[n2-1][n1-1] = 1.0;
    }
Therefore, the functions defined below do not allocate true
multi-dimensional arrays, as described in the C specification.
Instead, they allocate and initialize pointers (and pointers to
pointers) so that, for example, a[i2][i1] behaves like a 2-D array.

The array dimensions are numbered, which makes it easy to add
functions for arrays of higher dimensions.  In particular,
the 1st dimension of length n1 is always the fastest dimension,
the 2nd dimension of length n2 is the next fastest dimension,
and so on.  Note that the 1st (fastest) dimension n1 is the
first argument to the allocation functions defined below, but
that the 1st dimension is the last subscript in a[i2][i1].
(This is another important difference between C and Fortran.)

The allocation of pointers to pointers implies that more storage
is required than is necessary to hold a true multi-dimensional array.
The fraction of the total storage allocated that is used to hold
pointers is approximately 1/(n1+1).  This extra storage is unlikely
to represent a significant waste for large n1.

The functions defined below are significantly different from similar
functions described by Press et al, 1988, NR in C.
In particular, the functions defined below:
    (1) Allocate arrays of arbitrary size elements.
    (2) Allocate contiguous storage for arrays.
    (3) Return NULL if allocation fails (just like malloc).
    (4) Do not provide arbitrary lower and upper bounds for arrays.

Contiguous storage enables an allocated multi-dimensional array to
be passed to a C function that expects a one-dimensional array.
For example, to allocate and zero an n1 by n2 two-dimensional array
of floats, one could use
    a = alloc2(n1,n2,sizeof(float));
    zeroFloatArray(n1*n2,a[0]);
where zeroFloatArray is a function defined as
    void zeroFloatArray(int n, float *a)
    {
        int i;
        for (i=0; i<n; i++)
            a[i] = 0.0;
    }

Internal error handling and arbitrary array bounds, if desired,
should be implemented in functions that call the functions defined
below, with the understanding that these enhancements may limit
portability.

******************************************************************************
Author:    	Dave Hale, Colorado School of Mines, 12/31/89
                Zhaobo Meng, added 4D, 5D and 6D functions, 1996
*****************************************************************************/
/**************** end self doc ********************************/
/* allocate a 1-d array */
void *alloc1(size_t n1, size_t size)
{
    void *p;

    if ((p = malloc(n1 * size)) == NULL)
        err("%s : alloc1: malloc failed", __FILE__);
    return p;
}
/* re-allocate a 1-d array */
void *realloc1(void *v, size_t n1, size_t size)
{
    void *p;

    if ((p = realloc(v, n1 * size)) == NULL)
        err("%s : realloc1: realloc failed", __FILE__);
    return p;
}
/* free a 1-d array */
void free1(void *p)
{
    if (p)
    {
        free(p);
    }
}
/* allocate a 2-d array */
void **alloc2(size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p = (void **)malloc(n2 * sizeof(void *))) == NULL)
        err("%s : alloc2: malloc failed", __FILE__);
    if ((p[0] = (void *)malloc(n2 * n1 * size)) == NULL)
    {
        free(p);
        err("%s : alloc2: malloc failed", __FILE__);
    }
    for (i2 = 0; i2 < n2; i2++)
        p[i2] = (char *)p[0] + size * n1 * i2;
    return p;
}

/* free a 2-d array */
void free2(void **p)
{
    free(p[0]);
    free(p);
}
/* allocate a 3-d array */
void ***alloc3(size_t n1, size_t n2, size_t n3, size_t size)
{
    size_t i3, i2;
    void ***p;

    if ((p = (void ***)malloc(n3 * sizeof(void **))) == NULL)
        err("%s : alloc3: malloc failed", __FILE__);
    if ((p[0] = (void **)malloc(n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p);
        err("%s : alloc3: malloc failed", __FILE__);
    }
    if ((p[0][0] = (void *)malloc(n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0]);
        free(p);
        err("%s : alloc3: malloc failed", __FILE__);
    }

    for (i3 = 0; i3 < n3; i3++)
    {
        p[i3] = p[0] + n2 * i3;
        for (i2 = 0; i2 < n2; i2++)
            p[i3][i2] = (char *)p[0][0] + size * n1 * (i2 + n2 * i3);
    }
    return p;
}

/* free a 3-d array */
void free3(void ***p)
{
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 4-d array */
void ****alloc4(size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
    size_t i4, i3, i2;
    void ****p;

    if ((p = (void ****)malloc(n4 * sizeof(void ***))) == NULL)
        err("%s : alloc4: malloc failed", __FILE__);
    if ((p[0] = (void ***)malloc(n4 * n3 * sizeof(void **))) == NULL)
    {
        free(p);
        err("%s : alloc4: malloc failed", __FILE__);
    }
    if ((p[0][0] = (void **)malloc(n4 * n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p[0]);
        free(p);
        err("%s : alloc4: malloc failed", __FILE__);
    }
    if ((p[0][0][0] = (void *)malloc(n4 * n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0][0]);
        free(p[0]);
        free(p);
        err("%s : alloc4: malloc failed", __FILE__);
    }
    for (i4 = 0; i4 < n4; i4++)
    {
        p[i4] = p[0] + i4 * n3;
        for (i3 = 0; i3 < n3; i3++)
        {
            p[i4][i3] = p[0][0] + n2 * (i3 + n3 * i4);
            for (i2 = 0; i2 < n2; i2++)
                p[i4][i3][i2] = (char *)p[0][0][0] +
                                size * n1 * (i2 + n2 * (i3 + n3 * i4));
        }
    }
    return p;
}

/* free a 4-d array */
void free4(void ****p)
{
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 5-d array */
void *****alloc5(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
    size_t i5, i4, i3, i2;
    void *****p;

    if ((p = (void *****)malloc(n5 * sizeof(void ****))) == NULL)
        err("%s : alloc5: malloc failed", __FILE__);
    if ((p[0] = (void ****)malloc(n5 * n4 * sizeof(void ***))) == NULL)
    {
        free(p);
        err("%s : alloc5: malloc failed", __FILE__);
    }
    if ((p[0][0] = (void ***)malloc(n5 * n4 * n3 * sizeof(void **))) == NULL)
    {
        free(p[0]);
        free(p);
        err("%s : alloc5: malloc failed", __FILE__);
    }
    if ((p[0][0][0] = (void **)malloc(n5 * n4 * n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p[0][0]);
        free(p[0]);
        free(p);
        err("%s : alloc5: malloc failed", __FILE__);
    }
    if ((p[0][0][0][0] = (void *)malloc(n5 * n4 * n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        err("%s : alloc5: malloc failed", __FILE__);
    }
    for (i5 = 0; i5 < n5; i5++)
    {
        p[i5] = p[0] + n4 * i5;
        for (i4 = 0; i4 < n4; i4++)
        {
            p[i5][i4] = p[0][0] + n3 * (i4 + n4 * i5);
            for (i3 = 0; i3 < n3; i3++)
            {
                p[i5][i4][i3] = p[0][0][0] + n2 * (i3 + n3 * (i4 + n4 * i5));
                for (i2 = 0; i2 < n2; i2++)
                    p[i5][i4][i3][i2] = (char *)p[0][0][0][0] +
                                        size * n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * i5)));
            }
        }
    }
    return p;
}

/* free a 5-d array */
void free5(void *****p)
{
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 6-d array */
void ******alloc6(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                  size_t size)
{
    size_t i6, i5, i4, i3, i2;
    void ******p;

    if ((p = (void ******)malloc(n6 * sizeof(void *****))) == NULL)
        err("%s : alloc6: malloc failed", __FILE__);

    if ((p[0] = (void *****)malloc(n6 * n5 * sizeof(void ****))) == NULL)
    {
        free(p);
        err("%s : alloc6: malloc failed", __FILE__);
    }

    if ((p[0][0] = (void ****)malloc(n6 * n5 * n4 * sizeof(void ***))) == NULL)
    {
        free(p[0]);
        free(p);
        err("%s : alloc6: malloc failed", __FILE__);
    }
    if ((p[0][0][0] = (void ***)malloc(n6 * n5 * n4 * n3 * sizeof(void **))) == NULL)
    {
        free(p[0][0]);
        free(p[0]);
        free(p);
        err("%s : alloc6: malloc failed", __FILE__);
    }
    if ((p[0][0][0][0] = (void **)malloc(n6 * n5 * n4 * n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        err("%s : alloc6: malloc failed", __FILE__);
    }
    if ((p[0][0][0][0][0] = (void *)malloc(n6 * n5 * n4 * n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0][0][0][0]);
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        err("%s : alloc6: malloc failed", __FILE__);
    }

    for (i6 = 0; i6 < n6; i6++)
    {
        p[i6] = p[0] + n5 * i6;
        for (i5 = 0; i5 < n5; i5++)
        {
            p[i6][i5] = p[0][0] + n4 * (i5 + n5 * i6);
            for (i4 = 0; i4 < n4; i4++)
            {
                p[i6][i5][i4] = p[0][0][0] + n3 * (i4 + n4 * (i5 + n5 * i6));
                for (i3 = 0; i3 < n3; i3++)
                {
                    p[i6][i5][i4][i3] = p[0][0][0][0] + n2 * (i3 + n3 * (i4 + n4 * (i5 + n5 * i6)));
                    for (i2 = 0; i2 < n2; i2++)
                        p[i6][i5][i4][i3][i2] =
                            (char *)p[0][0][0][0][0] +
                            size * n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * (i5 + n5 * i6))));
                }
            }
        }
    }
    return p;
}

/* free a 6-d array */
void free6(void ******p)
{
    free(p[0][0][0][0][0]);
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 1-d array of ints */
int *alloc1int(size_t n1)
{
    return (int *)alloc1(n1, sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, size_t n1)
{
    return (int *)realloc1(v, n1, sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
    free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(size_t n1, size_t n2)
{
    return (int **)alloc2(n1, n2, sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(size_t n1, size_t n2, size_t n3)
{
    return (int ***)alloc3(n1, n2, n3, sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
    free3((void ***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(size_t n1)
{
    return (float *)alloc1(n1, sizeof(float));
}

/* re-allocate a 1-d array of floats */
float *realloc1float(float *v, size_t n1)
{
    return (float *)realloc1(v, n1, sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
    free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
    return (float **)alloc2(n1, n2, sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
    return (float ***)alloc3(n1, n2, n3, sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
    free3((void ***)p);
}

/* allocate a 4-d array of floats, added by Zhaobo Meng, 1997 */
float ****alloc4float(size_t n1, size_t n2, size_t n3, size_t n4)
{
    return (float ****)alloc4(n1, n2, n3, n4, sizeof(float));
}

/* free a 4-d array of floats, added by Zhaobo Meng, 1997 */
void free4float(float ****p)
{
    free4((void ****)p);
}

/* allocate a 5-d array of floats, added by Zhaobo Meng, 1997 */
float *****alloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
    return (float *****)alloc5(n1, n2, n3, n4, n5, sizeof(float));
}

/* free a 5-d array of floats, added by Zhaobo Meng, 1997 */
void free5float(float *****p)
{
    free5((void *****)p);
}

/* allocate a 6-d array of floats, added by Zhaobo Meng, 1997 */
float ******alloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6)
{
    return (float ******)alloc6(n1, n2, n3, n4, n5, n6, sizeof(float));
}

/* free a 6-d array of floats, added by Zhaobo Meng, 1997 */
void free6float(float ******p)
{
    free6((void ******)p);
}

/* allocate a 4-d array of ints, added by Zhaobo Meng, 1997 */
int ****alloc4int(size_t n1, size_t n2, size_t n3, size_t n4)
{
    return (int ****)alloc4(n1, n2, n3, n4, sizeof(int));
}

/* free a 4-d array of ints, added by Zhaobo Meng, 1997 */
void free4int(int ****p)
{
    free4((void ****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
int *****alloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
    return (int *****)alloc5(n1, n2, n3, n4, n5, sizeof(int));
}

/* free a 5-d array of ints, added by Zhaobo Meng, 1997 */
void free5int(int *****p)
{
    free5((void *****)p);
}

/* allocate a 5-d array of chars, added by Zhaobo Meng, 1997 */
unsigned char *****alloc5uchar(size_t n1, size_t n2, size_t n3,
                               size_t n4, size_t n5)
{
    return (unsigned char *****)alloc5(n1, n2, n3, n4, n5, sizeof(unsigned char));
}

/* free a 5-d array of chars, added by Zhaobo Meng, 1997 */
void free5uchar(unsigned char *****p)
{
    free5((void *****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
unsigned short *****alloc5ushort(size_t n1, size_t n2, size_t n3,
                                 size_t n4, size_t n5)
{
    return (unsigned short *****)alloc5(n1, n2, n3, n4, n5, sizeof(unsigned short));
}

/* free a 5-d array of shorts, added by Zhaobo Meng, 1997 */
void free5ushort(unsigned short *****p)
{
    free5((void *****)p);
}

/* allocate a 6-d array of ints, added by Zhaobo Meng, 1997 */
unsigned short ******alloc6ushort(size_t n1, size_t n2, size_t n3,
                                  size_t n4, size_t n5, size_t n6)
{
    return (unsigned short ******)alloc6(n1, n2, n3, n4, n5, n6, sizeof(unsigned short));
}

/* free a 6-d array of shorts, added by Zhaobo Meng, 1997 */
void free6ushort(unsigned short ******p)
{
    free6((void ******)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(size_t n1)
{
    return (double *)alloc1(n1, sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, size_t n1)
{
    return (double *)realloc1(v, n1, sizeof(double));
}

/* free a 1-d array of doubles */
void free1double(double *p)
{
    free1(p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(size_t n1, size_t n2)
{
    return (double **)alloc2(n1, n2, sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(size_t n1, size_t n2, size_t n3)
{
    return (double ***)alloc3(n1, n2, n3, sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
    free3((void ***)p);
}

/* allocate a 1-d array of cpxs */
cpx *alloc1cpx(size_t n1)
{
    return (cpx *)alloc1(n1, sizeof(cpx));
}

/* re-allocate a 1-d array of cpxs */
cpx *realloc1cpx(cpx *v, size_t n1)
{
    return (cpx *)realloc1(v, n1, sizeof(cpx));
}

/* free a 1-d array of cpxs */
void free1cpx(cpx *p)
{
    free1(p);
}

/* allocate a 2-d array of cpxs */
cpx **alloc2cpx(size_t n1, size_t n2)
{
    return (cpx **)alloc2(n1, n2, sizeof(cpx));
}

/* free a 2-d array of cpxs */
void free2cpx(cpx **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of cpxs */
cpx ***alloc3cpx(size_t n1, size_t n2, size_t n3)
{
    return (cpx ***)alloc3(n1, n2, n3, sizeof(cpx));
}

/* free a 3-d array of cpxs */
void free3cpx(cpx ***p)
{
    free3((void ***)p);
}

/*my code starts here to alloc double cpx functions*/

/* allocate a 1-d array of cpxs */
zpx *alloc1zpx(size_t n1)
{
    return (zpx *)alloc1(n1, sizeof(zpx));
}

/* re-allocate a 1-d array of cpxs */
zpx *realloc1zpx(zpx *v, size_t n1)
{
    return (zpx *)realloc1(v, n1, sizeof(zpx));
}

/* free a 1-d array of cpxs */
void free1zpx(zpx *p)
{
    free1(p);
}

/* allocate a 2-d array of cpxs */
zpx **alloc2zpx(size_t n1, size_t n2)
{
    return (zpx **)alloc2(n1, n2, sizeof(zpx));
}

/* free a 2-d array of cpxs */
void free2zpx(zpx **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of cpxs */
zpx ***alloc3zpx(size_t n1, size_t n2, size_t n3)
{
    return (zpx ***)alloc3(n1, n2, n3, sizeof(zpx));
}

/* free a 3-d array of cpxs */
void free3zpx(zpx ***p)
{
    free3((void ***)p);
}
/* added by Guo Jiadong add,2025 */
char *alloc1char(size_t n1)
{
    return (char *)alloc1(n1, sizeof(char));
}
char *realloc1char(char *v, size_t n1)
{
    return (char *)realloc1(v, n1, sizeof(char));
}
void free1char(char *p)
{
    free1(p);
}
/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************
SUBCALLS - routines for system functions with error checking

efopen		fopen with error check
efreopen	freopen with error check
efdopen		fdopen with error check
epopen		popen with error check
efclose		fclose with error check
epclose		pclose with error check
efflush		fflush with error check
eremove		remove with error check
erename		rename with error check
efseek		fseek with error check
efseeko		fseeko with error check
eftell		ftell with error check
eftello		ftello with error check
etmpfile	tmpfile with error check
erewind		rewind (dummy sub)
emkstemp	mkstemp with error check
emalloc		malloc with error check
erealloc	realloc with error check
ecalloc		calloc with error check
efgetpos	fgetpos with error check
efsetpos	fsetpos with error check
efread		fread with error check
efwrite		fwrite with error check

****************************************************************************
Function Prototypes:
FILE *efopen(const char *file, const char *mode);
FILE *efreopen(const char *file, const char *mode, FILE *stream1);
FILE *efdopen(int fd, const char *mode);
FILE *epopen(char *command, char *type);
int efclose(FILE *stream);
int epclose(FILE *stream);
int efflush(FILE *stream);
int eremove(const char *file);
int erename(const char *oldfile, const char* newfile);
int efseek(FILE *stream, off_t offset, int origin);
void erewind(FILE *stream);
long eftell(FILE *stream);
off_t eftello(FILE *streem);
int efseek(FILE *stream, off_t offset, int origin);
int efgetpos(FILE *stream, fpos_t *position);
int efsetpos(FILE *stream, const fpos_t *position);
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream);
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream);

****************************************************************************
Returns: All abort with message on error

efopen returns a FILE pointer
efreopen returns a FILE pointer
efdopen returns a FILE pointer
efclose returns 0
efflush returns 0
eremove returns 0
erename returns 0
efseek returns 0
efseeko returns 0
erewind is void
eftell returns file position indicator
eftello returns file position indicator
efgetpos returns 0
efsetpos returns 0
efread returns number of items actually read
efwrite returns number of items actually written

****************************************************************************
Notes:
Getting less than the number of bytes asked for on a fread
is *not* a system subroutine error--usually it just means
end of file has been reached.  However, it *might* be an error
in some application.  Similarly coming up empty is not a system
error, but might be an application error.  It is left to the user
to trap these instances.  In particular, feof can be used to test
for end of file.

****************************************************************************
References:
Rochkind, "Advanced UNIX Programming"
Kernighan and Pike, "The UNIX Programming Environment"
Kernighan and Ritchie, "The C Programming Language"
Mark Williams Company, "ANSI C--A Lexical Guide"

****************************************************************************
Authors: SEP: Rick Ottolini, Ron, Jon Claerbout, Stew Levin
CWP: Shuki Ronen, Jack Cohen
***************************************************************************/
/**************** end self doc ********************************/

FILE *efopen(const char *file, const char *mode)
{
    FILE *stream;

    if (NULL == (stream = fopen(file, mode)))
        err("%s: efopen: fopen failed", __FILE__);
    return stream;
}
FILE *efreopen(const char *file, const char *mode, FILE *stream1)
{
    FILE *stream2;

    if (NULL == (stream2 = freopen(file, mode, stream1)))
        err("%s: efreopen: freopen failed", __FILE__);

    return stream2;
}
FILE *efdopen(int fd, const char *mode)
{
    FILE *stream;

    if (NULL == (stream = (FILE *)fdopen(fd, mode)))
        err("%s: efdopen: fdopen failed", __FILE__);

    return stream;
}
FILE *epopen(char *command, char *type)
{
    FILE *stream;

    if (NULL == (stream = (FILE *)popen(command, type)))
        err("%s: epopen: popen failed", __FILE__);

    return stream;
}
int efclose(FILE *stream)
{
    int status;

    if (EOF == (status = fclose(stream)))
        err("%s: efclose: fclose failed", __FILE__);

    return status;
}
int epclose(FILE *stream)
{
    int status;

    if (EOF == (status = pclose(stream)))
        err("%s: epclose: pclose failed", __FILE__);

    return status;
}
int efflush(FILE *stream)
{
    int status;

    if (EOF == (status = fflush(stream)))
        err("%s: efflush: fflush failed", __FILE__);

    return status;
}
int eremove(const char *file)
{
    int status;

    if ((status = remove(file)))
        err("%s: eremove: remove failed", __FILE__);

    return status;
}
int erename(const char *oldfile, const char *newfile)
{
    int status;

    if ((status = rename(oldfile, newfile)))
        err("%s: erename: rename failed", __FILE__);

    return status;
}
int efseek(FILE *stream, off_t offset, int origin)
{
    if (fseek(stream, offset, origin)) /* non-zero => error */
        err("%s: efseek: fseek failed", __FILE__);

    return 0;
}
void erewind(FILE *stream) /* dummy function */
{
    rewind(stream);
    return;
}
long eftell(FILE *stream)
{
    long position;

    if (-1L == (position = ftell(stream)))
        err("%s: eftell: ftell failed", __FILE__);

    return position;
}
int efseeko(FILE *stream, off_t offset, int origin)
{

    /* non-zero => error */
    if (fseeko(stream, (off_t)offset, (int)origin))
        err("%s : efseeko: fseeko failed", __FILE__);

    return 0;
}
off_t eftello(FILE *streem)
{
    off_t eposition;
    off_t test = -1;

    eposition = ftello(streem);
    if (test == eposition)
    {
        fprintf(stderr, "sizeof(off_t)=%lu\n",
                (unsigned long)sizeof(eposition));
    }

    return eposition;
}
int efgetpos(FILE *stream, fpos_t *position)
{
    int status;

    if ((status = fgetpos(stream, position)))
        err("%s: efgetpos: fgetpos failed", __FILE__);

    return status;
}

int efsetpos(FILE *stream, const fpos_t *position)
{
    int status;

    if ((status = fsetpos(stream, position)))
        err("%s: efsetpos: fsetpos failed", __FILE__);

    return status;
}
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream)
{
    size_t nread;

    if (!size)
        err("%s: efread: fread given 0 item size", __FILE__);

    nread = fread(bufptr, size, count, stream);

    if (nread != count && ferror(stream))
        err("%s: efread: fread only %d items of %d",
            __FILE__, nread, count);

    return nread;
}
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream)
{
    size_t nwrite;

    nwrite = fwrite(bufptr, size, count, stream);

    if (nwrite != count)
        err("%s: efwrite: fwrite only %d items of %d",
            __FILE__, nwrite, count);

    return nwrite;
}

/*getpar.c  Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.			*/
/*********************** self documentation **********************/
/*****************************************************************************
GETPARS - Functions to GET PARameterS from the command line. Numeric
    parameters may be single values or arrays of int, uint,
    short, ushort, long, ulong, float, or double.  Single character
    strings (type string or char *) may also be gotten.
    Arrays of strings, delimited by, but not containing
    commas are permitted.

The functions are:

initargs 	Makes command line args available to subroutines (re-entrant).
        Every par program starts with this call!
getparint		get integers
getparuint		get unsigned integers
getparshort		get short integers
getparushort		get unsigned short integers
getparlong		get long integers
getparulong		get unsigned long integers
getparfloat		get float
getpardouble		get double
getparstring		get a single string
getparstringarray	get string array (fields delimited by commas)
getpar			get parameter by type
getnparint		get n'th occurrence of integer
getnparuint		get n'th occurrence of unsigned int
getnparshort		get n'th occurrence of short integer
getnparushort		get n'th occurrence of unsigned short int
getnparlong		get n'th occurrence of long integer
getnparulong		get n'th occurrence of unsigned long int
getnparfloat		get n'th occurrence of float
getnpardouble		get n'th occurrence of double
getnparstring		get n'th occurrence of string
getnparstringarray	get n'th occurrence of string array
getnpar			get n'th occurrence by type
countparname		return the number of times a parameter names is used
countparval		return the number of values in the last occurrence
                of a parameter
countnparval		return the number of values in the n'th occurrence
                of a parameter
getPar			Promax compatible version of getpar
checkpars()		check the argument list for typos
******************************************************************************
Function Prototypes:
void initargs (int argc, char **argv);
int getparint (char *name, int *p);
int getparuint (char *name, unsigned int *p);
int getparshort (char *name, short *p);
int getparushort (char *name, unsigned short *p);
int getparlong (char *name, long *p);
int getparulong (char *name, unsigned long *p);
int getparfloat (char *name, float *p);
int getpardouble (char *name, double *p);
int getparstring (char *name, char **p);
int getparstringarray (char *name, char **p);
int getnparint (int n, char *name, int *p);
int getnparuint (int n, char *name, unsigned int *p);
int getnparshort (int n, char *name, short *p);
int getnparushort (int n, char *name, unsigned short *p);
int getnparlong (int n, char *name, long *p);
int getnparulong (int n, char *name, unsigned long *p);
int getnparfloat (int n, char *name, float *p);
int getnpardouble (int n, char *name, double *p);
int getnparstring (int n, char *name, char **p);
int getnparstringarray (int n, char *name, char **p);
int getnpar (int n, char *name, char *type, void *ptr);
int countparname (char *name);
int countparval (char *name);
int countnparval (int n, char *name);
void getPar(char *name, char *type, void *ptr);
void checkpars( void );

******************************************************************************
Notes:
Here are some usage examples:

    ... if integer n not specified, then default to zero.
    if (!getparint("n", &n)) n = 0;

    ... if array of floats vx is specified, then
    if (nx=countparval("vx")) {
        ... allocate space for array
        vx = alloc1float(nx);
        ... and get the floats
        getparfloat("vx",vx);
    }

The command line for the above examples might look like:
    progname n=35 vx=3.21,4,9.5
    Every par program starts with this call!

More examples are provided in the DTEST code at the end of this file.

The functions: eatoh, eatou, eatol, eatov, eatoi, eatop used
below are versions of atoi that check for overflow.  The source
file for these functions is atopkge.c.

******************************************************************************
Authors:
Rob Clayton & Jon Claerbout, Stanford University, 1979-1985
Shuki Ronen & Jack Cohen, Colorado School of Mines, 1985-1990
Dave Hale, Colorado School of Mines, 05/29/90
Credit to John E. Anderson for re-entrant initargs 03/03/94
*****************************************************************************/
/**************** end self doc ********************************/

/* parameter table */
typedef struct
{
    char *name;     /* external name of parameter	*/
    char *asciival; /* ascii value of parameter	*/
} pointer_table;
/* global variables declared and used internally */
static pointer_table *argtbl; /* parameter table		*/
static int nargs;             /* number of args that parse	*/
static int tabled = false;    /* true when parameters tabled 	*/
static size_t targc;          /* total number of args		*/
static char **targv;          /* pointer to arg strings	*/
static char *argstr;          /* storage for command line	*/

/* functions declared and used internally */
static int getparindex(int n, char *name);
static void getparinit(void);
static void tabulate(size_t argc, char **argv);
static char *getpfname(void);
static int ccount(char c, char *s);
/*--------------------------------------------------------------------*\
  These variables are used by checkpars() to warn of parameter typos.
  par= does not use getpar() so we need to store that as the first
  parameter name.  lheaders= is buried in fgettr() so we intialize
  that also
  \*--------------------------------------------------------------------*/

#define PAR_NAMES_MAX 512

static char *par_names[PAR_NAMES_MAX];
static int par_count = 0;
static int parcheck = 0;
int xargc;
char **xargv;
/* make command line args available to subroutines -- re-entrant version */
void initargs(int argc, char **argv)
{
    /* added by Guo Jiadong,2025 */
    if (NULL != prog)
    {
        free(prog);
    }
    char *s = strrchr(argv[0], '/');
    s = (NULL != s) ? s + 1 : argv[0];
    int len = strlen(s) + 1;
    prog = alloc1char(len);
    strcpy(prog, s);

    /*------------------------------------------------------*/
    memset(par_names, 0, sizeof(par_names));
    par_names[0] = "par";
    par_names[1] = "lheader";
    par_count = 2;

    xargc = argc;
    xargv = argv;
    if (tabled == true)
    {
        free(argstr);
        free(targv);
        free(argtbl);
    }
    tabled = false;
    return;
}
void strchop(char *s, char *t)
/***********************************************************************
strchop - chop off the tail end of a string "s" after a "," returning
      the front part of "s" as "t".
      ************************************************************************
Notes:
Based on strcpy in Kernighan and Ritchie's C [ANSI C] book, p. 106.
************************************************************************
Author: CWP: John Stockwell and Jack K. Cohen, July 1995
***********************************************************************/
{

    while ((*s != ',') && (*s != '\0'))
    {
        *t++ = *s++;
    }
    *t = '\0';
}
/* functions to get values for the last occurrence of a parameter name */
int getparint(char *name, int *ptr)
{
    return getnpar(0, name, "i", ptr);
}
int getparuint(char *name, unsigned int *ptr)
{
    return getnpar(0, name, "p", ptr);
}
int getparshort(char *name, short *ptr)
{
    return getnpar(0, name, "h", ptr);
}
int getparushort(char *name, unsigned short *ptr)
{
    return getnpar(0, name, "u", ptr);
}
int getparlong(char *name, long *ptr)
{
    return getnpar(0, name, "l", ptr);
}
int getparulong(char *name, unsigned long *ptr)
{
    return getnpar(0, name, "v", ptr);
}
int getparfloat(char *name, float *ptr)
{
    return getnpar(0, name, "f", ptr);
}
int getpardouble(char *name, double *ptr)
{
    return getnpar(0, name, "d", ptr);
}
int getparstring(char *name, char **ptr)
{
    return getnpar(0, name, "s", ptr);
}
int getparstringarray(char *name, char **ptr)
{
    return getnpar(0, name, "a", ptr);
}
int getpar(char *name, char *type, void *ptr)
{
    return getnpar(0, name, type, ptr);
}

/* functions to get values for the n'th occurrence of a parameter name */
int getnparint(int n, char *name, int *ptr)
{
    return getnpar(n, name, "i", ptr);
}
int getnparuint(int n, char *name, unsigned int *ptr)
{
    return getnpar(n, name, "p", ptr);
}
int getnparshort(int n, char *name, short *ptr)
{
    return getnpar(n, name, "h", ptr);
}
int getnparushort(int n, char *name, unsigned short *ptr)
{
    return getnpar(n, name, "u", ptr);
}
int getnparlong(int n, char *name, long *ptr)
{
    return getnpar(n, name, "l", ptr);
}
int getnparulong(int n, char *name, unsigned long *ptr)
{
    return getnpar(n, name, "v", ptr);
}
int getnparfloat(int n, char *name, float *ptr)
{
    return getnpar(n, name, "f", ptr);
}
int getnpardouble(int n, char *name, double *ptr)
{
    return getnpar(n, name, "d", ptr);
}
int getnparstring(int n, char *name, char **ptr)
{
    return getnpar(n, name, "s", ptr);
}
int getnparstringarray(int n, char *name, char **ptr)
{
    return getnpar(n, name, "a", ptr);
}
int getnpar(int n, char *name, char *type, void *ptr)
{
    int i;      /* index of name in symbol table	*/
    int j;      /* index for par_names[]		*/
    int nval;   /* number of parameter values found	*/
    char *aval; /* ascii field of symbol		*/

    /*--------------------------------------------------------------------*\
      getpar gets called in loops reading traces in some programs.  So
      check for having seen this name before. Also make sure we don't
      walk off the end of the table.
      \*--------------------------------------------------------------------*/

    if (parcheck && strcmp("lheader", name))
    {
        fprintf(stderr, "getpar() call after checkpars(): %s\n", name);
    }

    for (j = 0; j < par_count; j++)
    {
        if (!strcmp(par_names[j], name))
        {
            break;
        }
    }

    if (j >= par_count && par_count < PAR_NAMES_MAX)
    {
        par_names[par_count++] = name;
    }

    if (par_count == PAR_NAMES_MAX)
    {
        fprintf(stderr, " %s exceeded PAR_NAMES_MAX %d \n", xargv[0], PAR_NAMES_MAX);
    }

    if (xargc == 1)
        return 0;
    if (!tabled)
        getparinit();         /* Tabulate command line and parfile */
    i = getparindex(n, name); /* Get parameter index */
    if (i < 0)
        return 0; /* Not there */

    if (0 == ptr)
    {
        fprintf(stderr, "%s: getnpar called with 0 pointer, type = %s\n", __FILE__, type);
    }

    /*
     * handle string type as a special case, since a string
     * may contain commas.
     */
    if (type[0] == 's')
    {
        *((char **)ptr) = argtbl[i].asciival;
        return 1;
    }

    /* convert vector of ascii values to numeric values */
    for (nval = 0, aval = argtbl[i].asciival; *aval; nval++)
    {
        switch (type[0])
        {
        case 'i':
            *(int *)ptr = eatoi(aval);
            ptr = (int *)ptr + 1;
            break;
        case 'p':
            *(unsigned int *)ptr = eatop(aval);
            ptr = (unsigned int *)ptr + 1;
            break;
        case 'h':
            *(short *)ptr = eatoh(aval);
            ptr = (short *)ptr + 1;
            break;
        case 'u':
            *(unsigned short *)ptr = eatou(aval);
            ptr = (unsigned short *)ptr + 1;
            break;
        case 'l':
            *(long *)ptr = eatol(aval);
            ptr = (long *)ptr + 1;
            break;
        case 'v':
            *(unsigned long *)ptr = eatov(aval);
            ptr = (unsigned long *)ptr + 1;
            break;
        case 'f':
            *(float *)ptr = eatof(aval);
            ptr = (float *)ptr + 1;
            break;
        case 'd':
            *(double *)ptr = eatod(aval);
            ptr = (double *)ptr + 1;
            break;
        case 'a':
        {
            char *tmpstr = "";
            tmpstr = alloc1(strlen(aval) + 1, 1);

            strchop(aval, tmpstr);
            *(char **)ptr = tmpstr;
            ptr = (char **)ptr + 1;
        }
        break;
        default:
            fprintf(stderr, "%s: invalid parameter type = %s",
                    __FILE__, type);
        }
        while (*aval++ != ',')
        {
            if (!*aval)
                break;
        }
    }
    return nval;
}

void checkpars(void)
{

    int i;
    int j;
    char buf[256];

    if (getparint("verbose", &i) && i == 1)
    {

#ifdef SUXDR
        fprintf(stderr, "Using Big Endian SU data format w/ XDR.\n");
#else
        fprintf(stderr, "Using native byte order SU data format w/o XDR.\n");
#endif
    }

    for (j = 1; j < xargc; j++)
    {

        for (i = 0; i < par_count; i++)
        {
            sprintf(buf, "%s=", par_names[i]);

            if (!strncmp(buf, xargv[j], strlen(buf)))
            {
                break;
            }
        }
        if (i == par_count && strchr(xargv[j], '='))
        {
            fprintf(stderr, "Unknown %s argument %s\n", xargv[0], xargv[j]);
        }
    }

    parcheck = 1;
}
/* return number of occurrences of parameter name */
int countparname(char *name)
{
    int i, nname;

    if (xargc == 1)
        return 0;
    if (!tabled)
        getparinit();
    for (i = 0, nname = 0; i < nargs; ++i)
        if (!strcmp(name, argtbl[i].name))
            ++nname;
    return nname;
}

/* return number of values in n'th occurrence of parameter name */
int countnparval(int n, char *name)
{
    int i;

    if (xargc == 1)
        return 0;
    if (!tabled)
        getparinit();
    i = getparindex(n, name);
    if (i >= 0)
        return ccount(',', argtbl[i].asciival) + 1;
    else
        return 0;
}

/* return number of values in last occurrence of parameter name */
int countparval(char *name)
{
    return countnparval(0, name);
}
/*
 * Return the index of the n'th occurrence of a parameter name,
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
static int getparindex(int n, char *name)
{
    int i;
    if (n == 0)
    {
        for (i = nargs - 1; i >= 0; --i)
            if (!strcmp(name, argtbl[i].name))
                break;
        return i;
    }
    else
    {
        for (i = 0; i < nargs; ++i)
            if (!strcmp(name, argtbl[i].name))
                if (--n == 0)
                    break;
        if (i < nargs)
            return i;
        else
            return -1;
    }
}

/* Initialize getpar */
static void getparinit(void)
{
    static char *pfname; /* name of parameter file		*/
    FILE *pffd = NULL;   /* file id of parameter file		*/
    int pflen;           /* length of parameter file in bytes	*/
    int parfile;         /* parfile existence flag		*/
    int argstrlen = 0;
    char *pargstr; /* storage for parameter file args	*/
    int nread = 0; /* bytes fread				*/
    int i, j;      /* counters				*/
    int start = true;
    int debug = false;
    int valencete = false;

    tabled = true; /* remember table is built		*/

    /* Check if xargc was initiated */

    if (!xargc)
        fprintf(stderr, "%s: xargc=%d -- not initiated in main\n", __FILE__, xargc);
    /* Space needed for command lines */

    for (i = 1, argstrlen = 0; i < xargc; i++)
    {
        argstrlen += strlen(xargv[i]) + 1;
    }

    /* Get parfile name if there is one */

    if ((pfname = getpfname()))
    {
        parfile = true;
    }
    else
    {
        parfile = false;
    }

    if (parfile)
    {
        pffd = fopen(pfname, "r");

        /* Get the length */
        fseek(pffd, 0, SEEK_END);

        pflen = ftell(pffd);

        rewind(pffd);
        argstrlen += pflen;
    }
    else
    {
        pflen = 0;
    }

    /*--------------------------------------------------------------------*\
      Allocate space for command line and parameter file. The pointer
      table could be as large as the string buffer, but no larger.

      The parser logic has been completely rewritten to prevent bad
      input from crashing the program.

      Reginald H. Beardsley			    rhb@acm.org
      \*--------------------------------------------------------------------*/

    argstr = (char *)alloc1(argstrlen + 1, 1);
    targv = (char **)alloc1((argstrlen + 1) / 4, sizeof(char *));

    if (parfile)
    {
        /* Read the parfile */

        nread = fread(argstr, 1, pflen, pffd);
        if (nread != pflen)
        {
            fprintf(stderr, "%s: fread only %d bytes out of %d from %s\n",
                    __FILE__, nread, pflen, pfname);
        }
        fclose(pffd);
    }

    /* force input to valid 7 bit ASCII */

    for (i = 0; i < nread; i++)
    {
        argstr[i] &= 0x7F;
    }

    /* tokenize the input */

    j = 0;

    for (i = 0; i < nread; i++)
    {

        /* look for start of token */

        if (start)
        {

            /* getpars.c:475: warning: subscript has type `char' */
            if (isgraph((int)argstr[i]))
            {
                targv[j] = &(argstr[i]);
                start = !start;
                j++;
            }
            else
            {
                argstr[i] = 0;
            }

            /* terminate token */

            /* getpars.c:487: warning: subscript has type `char' */
        }
        else if (!valencete && isspace((int)argstr[i]))
        {
            argstr[i] = 0;
            start = !start;
        }

        /* toggle valencete semaphore */

        if (argstr[i] == '\'' || argstr[i] == '\"')
        {
            valencete = !valencete;
        }
    }

    /* display all tokens */

    if (debug)
    {

        i = 0;
        while (i < j && targv[i] != 0)
        {
            if (strlen(targv[i]))
            {
                fprintf(stderr, "%d -> %s\n", i, targv[i]);
            }
            i++;
        }
    }

    /* discard non-parameter tokens */

    i = 0;
    targc = 0;
    while (i < j && targv[i] != 0)
    {
        if (strchr(targv[i], '='))
        {
            targv[targc] = targv[i];
            targc++;
        }
        i++;
    }

    /* Copy command line arguments */

    for (j = 1, pargstr = argstr + pflen + 1; j < xargc; j++)
    {
        strcpy(pargstr, xargv[j]);
        targv[targc++] = pargstr;
        pargstr += strlen(xargv[j]) + 1;
    }

    /* Allocate space for the pointer table */

    argtbl = (pointer_table *)alloc1(targc, sizeof(pointer_table));

    /* Tabulate targv */

    tabulate(targc, targv);

    return;
}
#define PFNAME "par="
/* Get name of parameter file */
static char *getpfname(void)
{
    int i;
    size_t pfnamelen;

    pfnamelen = strlen(PFNAME);
    for (i = xargc - 1; i > 0; i--)
    {
        if (!strncmp(PFNAME, xargv[i], pfnamelen) && strlen(xargv[i]) != pfnamelen)
        {
            return xargv[i] + pfnamelen;
        }
    }
    return NULL;
}

#define iswhite(c) ((c) == ' ' || (c) == '\t' || (c) == '\n')

/* Install symbol table */
static void tabulate(size_t argc, char **argv)
{
    int i;
    char *eqptr;
    int debug = false;

    for (i = 0, nargs = 0; i < argc; i++)
    {
        eqptr = strchr(argv[i], '=');
        if (eqptr)
        {
            argtbl[nargs].name = argv[i];
            argtbl[nargs].asciival = eqptr + 1;
            *eqptr = (char)0;

            /* Debugging dump */
            if (debug)
            {
                fprintf(stderr,
                        "argtbl[%d]: name=%s asciival=%s\n",
                        nargs, argtbl[nargs].name, argtbl[nargs].asciival);
            }
            nargs++;
        }
    }
    return;
}

/* Count characters in a string */
static int ccount(char c, char *s)
{
    int i, count;
    for (i = 0, count = 0; s[i] != 0; i++)
        if (s[i] == c)
            count++;
    return count;
}