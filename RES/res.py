import numpy as np
import matplotlib.pyplot as plt
def read_wfd(file,nx,nz):
    with open(file,"rb") as f:
        data=np.fromfile(f,dtype=np.float32)
        return data.reshape(nx,nz).T


plt.rcParams['font.sans-serif'] = ['Nimbus Roman']
plt.rcParams['axes.unicode_minus'] = False  # 设置支持负号显示
plt.rcParams['svg.fonttype'] = 'none'
    
nx=400
nz=400
dh=10
res_cmap='seismic'
data=read_wfd("PSM.bin",nx,nz)
data2=read_wfd("FDM2.bin",nx,nz)
data4=read_wfd("FDM4.bin",nx,nz)
data6=read_wfd("FDM6.bin",nx,nz)
data8=read_wfd("FDM8.bin",nx,nz)
data10=read_wfd("FDM10.bin",nx,nz)
res2=data-data2
res4=data-data4
res6=data-data6
res8=data-data8
res10=data-data10


# res2=np.abs(res2)
# res4=np.abs(res4)
# res6=np.abs(res6)
# res8=np.abs(res8)
# res10=np.abs(res10)
bnd=0.1
fig=plt.figure()
ax=fig.add_subplot(111)
ax.imshow(data,cmap='seismic',vmax=0.5,vmin=-0.5, interpolation='bicubic',extent=[0 ,(nx-1)*dh ,(nz-1)*dh ,0],origin='upper',)
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.savefig("PSM.svg",dpi=300)
########################################################################
fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(res2,cmap=res_cmap,vmax=bnd,vmin=-bnd, interpolation='bicubic',extent=[0 ,(nx-1)*dh ,(nz-1)*dh ,0],origin='upper',)
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.colorbar()
plt.savefig("FDM2.svg",dpi=300)
########################################################################
fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(res4,cmap=res_cmap,vmax=bnd,vmin=-bnd, interpolation='bicubic',extent=[0 ,(nx-1)*dh ,(nz-1)*dh ,0],origin='upper',)
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.colorbar()
plt.savefig("FDM4.svg",dpi=300)
########################################################################
fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(res6,cmap=res_cmap,vmax=bnd,vmin=-bnd, interpolation='bicubic',extent=[0 ,(nx-1)*dh ,(nz-1)*dh ,0],origin='upper',)
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.colorbar()
plt.savefig("FDM6.svg",dpi=300)
########################################################################
fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(res8,cmap=res_cmap,vmax=bnd,vmin=-bnd, interpolation='bicubic',extent=[0 ,(nx-1)*dh ,(nz-1)*dh ,0],origin='upper',)
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.colorbar()
plt.savefig("FDM8.svg",dpi=300)
########################################################################
fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(res10,cmap=res_cmap,vmax=bnd,vmin=-bnd, interpolation='bicubic',extent=[0 ,(nx-1)*dh ,(nz-1)*dh ,0],origin='upper',)
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.colorbar()
plt.savefig("FDM10.svg",dpi=300)
########################################################################
plt.show()







