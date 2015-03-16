import h5py 
import sys
import pylab as P

datafile = sys.argv[-1]
data = h5py.File(datafile,"r")

img_extent = [data['/scales/sim_time'][0],data['/scales/sim_time'][-1],data['/scales/y/1'][0],data['/scales/y/1'][-1]]
fig = P.figure(figsize=(9,9))

# <ux**2>_x 
ux_ax = fig.add_axes([0.1,0.05,0.7,0.3])
ux_cbax = fig.add_axes([0.8,0.05,0.03,0.3])
ux_ax.set_xlabel('time',size=18)
ux = data['/tasks/<x kin en density>_x'][:,0,]
ux_im = ux_ax.imshow(ux.T,cmap='copper',extent=img_extent,aspect='auto',vmax=ux.max()*0.5)
uxcb = fig.colorbar(ux_im, cax=ux_cbax)
uxcb.set_label(r'$<u_x^2>_x$')

# <uy**2>_x 
uy_ax = fig.add_axes([0.1,0.35,0.7,0.3])
uy_cbax = fig.add_axes([0.8,0.35,0.03,0.3])
uy_ax.xaxis.set_visible(False)
uy = data['/tasks/<y kin en density>_x'][:,0,]
uy_im = uy_ax.imshow(uy.T,cmap='copper',extent=img_extent,aspect='auto')
fig.colorbar(uy_im, cax=uy_cbax)
uy_cbax.set_ylabel(r'$<u_y^2>_x$')

# vorticity
vort_ax = fig.add_axes([0.1,0.65,0.7,0.3])
vort_cbax = fig.add_axes([0.8,0.65,0.03,0.3])
vort_ax.xaxis.set_visible(False)
vort = data['/tasks/<vorticity>_x'][:,0,]
vort_im = vort_ax.imshow(vort.T,'RdBu',extent=img_extent,aspect='auto')
fig.colorbar(vort_im, cax=vort_cbax)
vort_cbax.set_ylabel(r'$<\omega>_x$')

fig.savefig('x_averages.png')
