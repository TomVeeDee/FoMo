from readfomo import readgoftcube, regulargoftcube, gaussfitgoftcube
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

if 'txtdata' not in locals(): txtdata=readgoftcube('../fomo-output.txt.gz')
if 'datdata' not in locals(): datdata=readgoftcube('../fomo-output.dat.gz')

emiss,xvec,yvec,lvec=regulargoftcube(txtdata)
# pl.imshow(emiss[50,:,:],extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),aspect='auto')
# pl.contourf(xvec,yvec,emiss[50,:,:],100,aspect=2./15.)
# pl.show()

if 'peak' not in locals(): peak,doppler,sigma,chisq,intens=gaussfitgoftcube(emiss,lvec)


pl.imshow(peak,extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),cmap=cm.hot,aspect='auto')
pl.title('peak intensity')
cb=pl.colorbar()
cb.set_label(r'$ergs\ cm^{-2} s^{-1} sr^{-1} \AA{}^{-1}$')
pl.show()
pl.imshow(intens,extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),cmap=cm.hot,aspect='auto')
pl.title('integrated intensity')
cb=pl.colorbar()
cb.set_label(r'$ergs\ cm^{-2} s^{-1} sr^{-1}$')
pl.show()
pl.imshow(doppler,extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),cmap='bwr',aspect='auto')
pl.title('Doppler shift')
cb=pl.colorbar()
cb.set_label('km/s')
pl.show()
pl.imshow(sigma,extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),cmap='brg',aspect='auto')
pl.title('spectral line width')
cb=pl.colorbar()
cb.set_label(r'$\sigma_w (\AA{})$')
pl.show()
pl.imshow(chisq,extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),cmap=cm.hot,aspect='auto',norm=LogNorm())
pl.title(r'$\chi^2$'+'(deviation from Gaussian)')
pl.colorbar()
pl.show()


# showing the peak intensity in a 3D way, like IDL surface plot
x,y=np.meshgrid(xvec,yvec)
fig = pl.figure(figsize=(18,8))
ax = fig.gca(projection='3d')
ax.view_init(elev=30.,azim=60)
surf=ax.plot_surface(x,y,peak,cmap=cm.hot,rstride=1,cstride=1)