import yt
import matplotlib.pyplot as plt

for i in range(10):
    ds = yt.load('./data/plt%05d' %i)
    ad0 = ds.covering_grid(level=0,
                           left_edge=ds.domain_left_edge,
                           dims=ds.domain_dimensions)
    Ex = ad0['Ex'].to_ndarray()
    plt.imshow( Ex[int(len(Ex)/2)], vmin=-0.5, vmax=0.5 )
    plt.savefig('data/img%05d.png' %i)
