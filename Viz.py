import yt
import matplotlib.pyplot as plt

import os, shutil

if os.path.exists('png'):
    shutil.rmtree('png')
os.mkdir('png')

for i in range(20):
    ds = yt.load('./data/plt%05d' %i)
    ad0 = ds.covering_grid(level=0,
                           left_edge=ds.domain_left_edge,
                           dims=ds.domain_dimensions)
    Ex = ad0['Ex'].to_ndarray()
    Nx, Ny, Nz = Ex.shape

    plt.figure(figsize=(10,5))
    plt.subplot(131)
    plt.imshow( Ex[int(Nx/2),:,:], vmin=-0.01, vmax=0.01 )
    plt.title('y-z plane')
    plt.subplot(132)
    plt.imshow( Ex[:,int(Ny/2),:], vmin=-0.01, vmax=0.01 )
    plt.title('x-z plane')
    plt.subplot(133)
    plt.imshow( Ex[:,:,int(Nz/2)], vmin=-0.01, vmax=0.01 )
    plt.title('x-y plane')
    
    plt.savefig('png/img%05d.png' %i)
