# High-precision acoustic modeling with second-order staggered difference

![](https://img.shields.io/badge/License-GPLv3-blue)  ![](https://img.shields.io/badge/Author-Jiadong_Guo-blue)  ![](https://img.shields.io/badge/Email-jdongguo@126.com-blue)  ![](https://img.shields.io/badge/Language-C_Shell_Python-blue)  ![](https://img.shields.io/badge/System-Linux-blue)  ![](https://img.shields.io/badge/Dependencies-MPI_OpenBlas-blue)

This program realizes second-order interleaved differential acoustic wave simulation, using PML boundary conditions, comparing the role of regular mesh and interleaved mesh.

## References

- Du Z, Liu J, Liu J, Xu F, Li Y. High-precision acoustic modeling with second-order staggered difference[J]. Arabian Journal of Geosciences, 2017, 10(21): 473.

## Code structure

- cstd.c: common function library

- ricker.c: ricker wavelet

- coeff.c: calculated difference coefficient

- regular.c: regular grid simulation

- sg.c: staggered grid simulation


## Instructions to run

1. run scons: compile
1. run regmodeling: regular grid simulation
1. run sgmodeling: staggered grid simulation
