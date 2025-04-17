import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# 读取二进制文件
def read_binary_float_file(filename, nx, nz):
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)  # 读取float32数据
    return data.reshape((nx, nz)).T

# 主程序
if len(sys.argv)<7:
    print("使用方法: python plot.py 波场文件名 nx nz dx dz vmax")
    sys.exit(1)
filename = sys.argv[1] # 替换为你的文件名
nx=int(sys.argv[2])
nz=int(sys.argv[3])
dx=float(sys.argv[4])
dz=float(sys.argv[5])
vmax=float(sys.argv[6])
data = read_binary_float_file(filename,nx,nz)

# 创建图形
fig=plt.figure(figsize=(8, 6))
ax=fig.add_subplot(1,1,1)
ax.imshow(data, 
           cmap='seismic',
           interpolation='bicubic',
           aspect='equal',
           origin='lower',
           vmax=vmax,vmin=-vmax,
           extent=[0 ,(nx-1)*dx, 0 ,(nz-1)*dz]
           )
ax.set_title('Wavefield')
ax.set_xlabel('Distance(m)')
ax.set_ylabel('Depth(m)')
ax.axis([0,(nx-1)*dx,(nz-1)*dz,0])
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('top')
plt.show()
