import matplotlib.animation as animation
import numpy as np
import pylab as plt
import matplotlib.cm as cm
from readfomo import readgoftcube, readgoftcube_units_chianti, regulargoftcube, gaussfitgoftcube
import glob
import sys
import os
import sunpy.visualization.colormaps.cm

# store filenames as a sorted list
unsorted=glob.glob(sys.argv[1])
filelist=sorted(unsorted)

# create data cube from file list
print("Reading file:",filelist[0])
data,units,chiantifile=readgoftcube_units_chianti(filelist[0]) #modified by Krishna
emiss,xvec,yvec,lvec=regulargoftcube(data)
allemiss=emiss
for filename in filelist[1:]:
    print("Reading file:",filename)
    data=readgoftcube(filename)
    emiss,xvec,yvec,lvec=regulargoftcube(data)
    allemiss=np.dstack((allemiss,emiss))

# initiate plot
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
ax.set_xlabel('x ('+units[0]+')')
ax.set_ylabel('y ('+units[1]+')')

# set colour map
cmap=cm.hot
aiafilter=os.path.basename(chiantifile)
aiafilter=aiafilter.decode('ascii') #added by Krishna
if ('goft_table_aia094' in aiafilter): cmap=plt.get_cmap('sdoaia094')
if ('goft_table_aia131' in aiafilter): cmap=plt.get_cmap('sdoaia131')
if ('goft_table_aia171' in aiafilter): cmap=plt.get_cmap('sdoaia171')
if ('goft_table_aia193' in aiafilter): cmap=plt.get_cmap('sdoaia193')
if ('goft_table_aia211' in aiafilter): cmap=plt.get_cmap('sdoaia211')
if ('goft_table_aia304' in aiafilter): cmap=plt.get_cmap('sdoaia304')
if ('goft_table_aia335' in aiafilter): cmap=plt.get_cmap('sdoaia335')

# cut off the emission at percent_emiss of maximum intensity of all the datacubes
# even if it is 99%, it clips off very bright pixels, allowing for the entire field of view to be rendered
percent_emiss=99
im = ax.imshow(allemiss[:,:,0],extent=(np.amin(xvec),np.amax(xvec),np.amin(yvec),np.amax(yvec)),cmap=cmap,aspect='auto',vmax=np.percentile(allemiss,percent_emiss),vmin=np.amin(allemiss))
annotation=ax.annotate(format(0,'04d'),(.93,.97),textcoords='axes fraction',color='white')
plt.title('integrated intensity')
cb=plt.colorbar(im)
cb.set_label('DN/s')

plt.tight_layout()

def init():
    return im,annotation

def update_img(filename):
    im.set_data(allemiss[:,:,filename])
    annotation.set_text(format(filename,'04d'))
    return im,annotation

# animate the frames
dpi=300
ani = animation.FuncAnimation(fig,update_img,frames=len(filelist),interval=30,blit=False,init_func=init)
writer = animation.writers['ffmpeg'](fps=5)
ani.save('demo.mp4',writer=writer,dpi=dpi)
    
