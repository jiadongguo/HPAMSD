/*
求解一阶交错差分系数
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
        for (int i = 1; i <= n; i++)
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
#ifdef __TEST__
int main(int argc, char **argv)
{
    int n;
    scanf("%d", &n);
    printf("n=%d\n", n);
    float *a = (float *)malloc(sizeof(float) * n);
    printf("start\n");
    get_sg_coeff(n, a);
    printf("end\n");
    for (int i = 0; i < n; i++)
    {
        printf("%e\n", a[i]);
    }

    return 0;
}
#endif