nz=400 
nx=400
dh=10
fm=20
nt=200
dt=0.001
t0=0.1
c=3000

all:PSM.bin FDM2.bin FDM4.bin FDM6.bin FDM8.bin FDM10.bin 
PSM.bin:PSM
	./PSM nz=$(nz) nx=$(nx) dh=$(dh) fm=$(fm) nt=$(nt) dt=$(dt) t0=$(t0) c=$(c) wfd=$@
FDM2.bin:FDM
	./FDM nz=$(nz) nx=$(nx) dh=$(dh) fm=$(fm) nt=$(nt) dt=$(dt) t0=$(t0) c=$(c) order=2 wfd=$@
FDM4.bin:FDM
	./FDM nz=$(nz) nx=$(nx) dh=$(dh) fm=$(fm) nt=$(nt) dt=$(dt) t0=$(t0) c=$(c) order=4 wfd=$@
FDM6.bin:FDM
	./FDM nz=$(nz) nx=$(nx) dh=$(dh) fm=$(fm) nt=$(nt) dt=$(dt) t0=$(t0) c=$(c) order=6 wfd=$@
FDM8.bin:FDM
	./FDM nz=$(nz) nx=$(nx) dh=$(dh) fm=$(fm) nt=$(nt) dt=$(dt) t0=$(t0) c=$(c) order=8 wfd=$@
FDM10.bin:FDM
	./FDM nz=$(nz) nx=$(nx) dh=$(dh) fm=$(fm) nt=$(nt) dt=$(dt) t0=$(t0) c=$(c) order=10 wfd=$@





PSM:cfft2.c cfft2.h cstd.c PSM.c ricker.c
	gcc $^ -o $@ -lm -lfftw3f
FDM:cfft2.c cfft2.h cstd.c FDM.c ricker.c
	gcc $^ -o $@ -lm -lfftw3f
