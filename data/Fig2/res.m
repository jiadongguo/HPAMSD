clc
clear
close all
fd=fopen('1a.bin');
data=fread(fd,[400 400],'float');
fclose(fd);


fd=fopen('wfd1.bin');
data1=fread(fd,[400 400],'float');
fclose(fd);

r=data-data1;
rel=r./(data+1e-7);
max(max(abs(rel)))

imagesc(rel*100);
clim([-1 1]);
colorbar;
colormap('jet');