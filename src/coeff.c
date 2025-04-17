#include "cstd.h"
void get_sg_coeff(int n, float *a)
{
    float tmp1, tmp2;
    float wgt;
    int i, j;
    for (int m = 0; m < n; m++)
    {
        j = m + 1;
        tmp1 = 0;
        tmp2 = 0;
        for (i = 1; i <= n; i++)
        {
            if (i == j)
                continue;
            tmp1 += 2 * log(2 * i - 1);
            tmp2 += log(fabs((2 * i - 1) * (2 * i - 1) - (2 * j - 1) * (2 * j - 1)));
        }
        if (j % 2 == 0)
            wgt = -1;
        else
            wgt = 1;
        wgt /= 2 * j - 1;
        a[m] = exp(tmp1 - tmp2) * wgt;
    }
}
/*
求解一阶规则差分系数
 */
void get_reg_coeff(int n, float *a)
{
    float tmp1, tmp2;
    float wgt;
    int i, j;
    for (int m = 0; m < n; m++)
    {
        j = m + 1;
        tmp1 = 0;
        tmp2 = 0;
        for (i = 1; i <= n; i++)
        {
            if (i == j)
                continue;
            tmp1 += 2 * log(i);
            tmp2 += log(fabs(i * i - j * j));
        }

        if (j % 2 == 0)
            wgt = -1;
        else
            wgt = 1;
        wgt /= 2 * j;
        a[m] = wgt * exp(tmp1 - tmp2);
    }
}
void get_center_coeff2(int n, float *a)
{
    float tmp1, tmp2;
    float wgt;
    int i, j;
    for (int m = 0; m < n; m++)
    {
        j = m + 1;
        tmp1 = 0;
        tmp2 = 0;
        for (i = 1; i <= n; i++)
        {
            if (i == j)
                continue;
            tmp1 += 2 * log(i);
            tmp2 += log(fabs(i * i - j * j));
        }

        if (j % 2 == 0)
            wgt = -1;
        else
            wgt = 1;
        wgt /= j * j;
        a[m] = wgt * exp(tmp1 - tmp2);
    }
}
void get_sg_coeff2(int n, float *a)
{
    memset(a, 0, sizeof(float) * 2 * n);
    float *c1 = (float *)malloc(sizeof(float) * n);
    get_sg_coeff(n, c1);
    for (int i = 0; i < n; i++)
    {
        a[0] -= SQUARE(c1[i]);
    }
    for (int i = 1; i < n; i++)
    {
        for (int k = 1; k <= n; k++)
        {
            a[i] += c1[k - 1] * c1[i - k];
            a[i + n] += c1[k - 1] * c1[i - k];
        }
        for (int k = 1; k <= n - i; k++)
        {
            a[i] -= 2 * c1[k - 1] * c1[i + k - 1];
        }
    }
    free(c1);
}
#ifdef __TEST__
int main(int argc, char **argv)
{
    int n;
    scanf("%d", &n);
    printf("n=%d\n", n);
    float *a = (float *)malloc(sizeof(float) * n);
    float *b = (float *)malloc(sizeof(float) * n);
    printf("start\n");
    get_reg_coeff(n, a);
    get_sg_coeff(n, a);
    printf("end\n");
    for (int i = 0; i < n; i++)
    {
        printf("%e\n", a[i]);
        printf("%e\n", b[i]);
    }
    return 0;
}
#endif