import yt
import matplotlib.pyplot as plt

import os, shutil

if os.path.exists('png'):
    shutil.rmtree('png')
os.mkdir('png')


def plot_field(field, line, nlines, ad0, vmax):
    F = ad0[field].to_ndarray()
    Nx, Ny, Nz = F.shape
    if vmax is not None:
        vmin = -vmax
    else:
        vmin = None
    plt.subplot2grid( (nlines,3), (line,0))
    plt.imshow( F[int(Nx/2),:,:], vmin=vmin, vmax=vmax, origin='lower' )
#    plt.colorbar()
    plt.title(field + ': y-z plane')
    plt.subplot2grid( (nlines,3), (line,1))
    plt.imshow( F[:,int(Ny/2),:], vmin=vmin, vmax=vmax, origin='lower' )
#    plt.colorbar()
    plt.title(field + ': x-z plane')
    plt.subplot2grid( (nlines,3), (line,2))
    plt.imshow( F[:,:,int(Nz/2)], vmin=vmin, vmax=vmax, origin='lower' )
#    plt.colorbar()
    plt.title(field + ': x-y plane')


for i in range(20):
    ds = yt.load('./data/plt%05d' %i)
    ad0 = ds.covering_grid(level=0,
                           left_edge=ds.domain_left_edge,
                           dims=ds.domain_dimensions)
    nlines = 4
    plt.figure(figsize=(10,nlines*4))
    plot_field('Ex', 0, nlines, ad0, vmax=0.1)
    plot_field('Ey', 1, nlines, ad0, vmax=0.4)
    plot_field('Ez', 2, nlines, ad0, vmax=0.4)
    plot_field('rho1', 3, nlines, ad0, vmax=None)
#    plot_field('rho2', 4, nlines, ad0, vmax=None)
    plt.savefig('png/img%05d.png' %i)
