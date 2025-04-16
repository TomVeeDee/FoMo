import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import matplotlib.cm as cm
import imageio.v2 as imageio
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mplcursors
import csv  
import cv2

custom_path = os.path.expanduser('~/Users/dariasorokina/Desktop/FoMo/my_script/')
sys.path.append(custom_path)
from readfomo import readgoftcube, readgoftcubechianti, regulargoftcube, gaussfitgoftcube
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm, PowerNorm 
from matplotlib.patches import Circle

# Directory containing the FoMo output data files
data_directory = "/Users/dariasorokina/Desktop/FoMo/my_script/eclipse_paper/500x500"
output_directory = "Figures_eclipse_earth_new"
os.makedirs(output_directory, exist_ok=True)

# List all .dat.gz files in the directory
file_list = glob.glob(os.path.join(data_directory, '*.dat.gz'))
# Sort file list to process files in numerical order
file_list.sort(key=lambda x: int(''.join(filter(str.isdigit, os.path.basename(x)))))
print(f"Found {len(file_list)} files in {data_directory}")

image_files = []

for filename in file_list:
    try:
        # Create data cube from file list
        print(f'Reading file {filename}')
        data = readgoftcube(filename)
        emiss, xvec, yvec, lvec = regulargoftcube(data)
        X, Y = np.meshgrid(xvec, yvec)
        radius = np.sqrt(X**2 + Y**2)        
 
        # Store emission data from the first file to subtract from others
        if filename == file_list[0]:
            emiss_reference = emiss  # Store emission data from the first file
            vec_ref=xvec
            yvec_ref=yvec
            print(f'Emiss variable for reference file {filename}: {emiss_reference}')
            
        print(f'Emiss variable for {filename}: {emiss}')
        # Replace infinities and NaNs with zero
        allemiss=np.nan_to_num(emiss)
        #-emiss_reference )
        
        allemiss[allemiss <= 1e-20] = 1.0
        #self-consistently with FoMo convert to the MSB units
        allemiss_MSB=allemiss*0.36e-10

        # Initiate plot
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        ax.set_aspect('auto')
        #ax.set_xlabel('x (Mm)', fontsize=24)
        #ax.set_ylabel('y (Mm)', fontsize=24)
        
        # Create a mask for the circular region
        X, Y = np.meshgrid(xvec, yvec)
        radius = 6*695.7
        radius_ins= 0.9*695.7
        circle = (X**2 + Y**2) <= radius**2
        circle_ins = (X**2 + Y**2) >= radius_ins**2 

        final_mask = circle & circle_ins
        masked_allemiss_MSB = np.ma.array(allemiss_MSB, mask=~final_mask)
        #*(np.sqrt(X**2 + Y**2) ** 2)

        #im = ax.imshow( masked_allemiss_MSB, extent=(np.amin(xvec), np.amax(xvec), np.amin(yvec), np.amax(yvec)), 
        #                aspect='equal', norm=LogNorm(vmin=0.01, vmax=0.2), cmap='gray')
                     
        im = ax.imshow( masked_allemiss_MSB, extent=(np.amin(xvec), np.amax(xvec), np.amin(yvec), np.amax(yvec)), 
                        aspect='equal', norm=LogNorm(vmin= masked_allemiss_MSB.min(), vmax= masked_allemiss_MSB.max()), cmap='gray')
     
        plt.title('FoMo: Polarised Brightness')
        #Create a divider for the existing axes instance
        #divider = make_axes_locatable(ax)
        # Append axes to the right of ax, with 5% width of ax
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        # Create a colorbar with the new axes
        #cb = plt.colorbar(im, cax=cax)
        #cb.ax.tick_params(labelsize=18) #18
        #cb.set_label('MSB', fontsize=18)

        # Optional: Set axes limits and circle patches
        #ax.set_xlim(1*695.7, 14*695.7)
        #ax.set_ylim(-3*695.7, 3*695.7)
        ax.set_xlim(-5*695.7, 5*695.7)
        ax.set_ylim(-3.3*695.7, 3.3*695.7)
        plt.xticks([])  # Remove x-axis ticks
        plt.yticks([])  # Remove y-axis ticks
        plt.axis('off')  # Turn off the axes
        #ax.tick_params(axis='both', which='major', labelsize=14) #18
        circ1 = Circle((0, 0), 1.1*695.7, facecolor='k', edgecolor='k', lw=1)
        circ2 = Circle((0, 0), 6.0*695.7, facecolor='none', edgecolor='w', lw=1)
        
        ax.add_patch(circ1)
        
        cursor = mplcursors.cursor(im, hover=True)
        
        # Save the plot
        output_filename = os.path.join(output_directory, os.path.basename(filename).replace('.gz', '.png'))
        plt.savefig(output_filename, format='png')
        image_files.append(output_filename)
        print(f'Saved plot to {output_filename}')
        plt.show()
        plt.close(fig)

        # Add to image list for video
        image_files.append(output_filename)
        
    except Exception as e:
        print(f"Failed to process file {filename}: {e}")

# Sort image files list to ensure images are in the correct order for video creation
image_files.sort(key=lambda x: int(''.join(filter(str.isdigit, os.path.basename(x)))))
        
video_filename = os.path.join(output_directory, 'output_video.mp4')
with imageio.get_writer(video_filename, mode='I', macro_block_size=1) as writer:
    for file in image_files:
        image = imageio.imread(file)
        writer.append_data(image)

print(f'Created video: {video_filename}')
