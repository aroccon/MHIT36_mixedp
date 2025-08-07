# # %%
# %%html
# <style>
# .cell-output-ipywidget-background {
#    background-color: transparent !important;
# }
# .jp-OutputArea-output {
#    background-color: transparent;
# }  
# .dataframe th {
#     font-size: 6px;
# }
# .dataframe td {
#     font-size: 6px;
# }
# </style>

# %%
# %matplotlib widget
# %matplotlib inline

import ipympl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
# from ing_theme_matplotlib import mpl_style
# from qbstyles import mpl_style

import os
import glob
import cv2

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata


from jupyterthemes import jtplot

from matplotlib.widgets import Cursor

import h5py

jtplot.style(theme='monokai', context='notebook', ticks=True, grid=False)
plt.style.use("dark_background")
# mpl_style(dark=True)


#################################################################
foldername = '../multi/output/'

# PARAMETERS:

# fields to plot
fld = 'p'

# number of points in each direction
nx = 256

# compute the derivative of the fields (show them instead of the neormal fields)
# 0: no derivative
# 1: x derivative
# 2: y derivative
# 3: z derivative
# list more flag to compute consecutive derivatives (forder 1 FD)
derivative_vec = [0]

# # normal direction of the 2D slice:
# 1: x-direction
# 2: y-direction
# 3: z-direction
slice_dir = 1

# index to take the slice (from 1 to nx_i, choose -1 for computing the average)
slice_idx = 128

# slice_idx = 222

slice_idx = round(nx/2)

# time_steps to plot 
# ts_vec = [0]
ts_vec = range(800000,840500,10000)
ts_vec = range(0,151,10)

# ts_vec = [5]

# value for the fontsize:
fontsize_val = 10

# slice of slice (leave -1 to compute the mean)
# meanprof_slice = -1

meanprof_slice = round(nx/2)


# value of the figure size
figsize_val = 10

# save figure in png format
savefig_flag = 0

# Set your desired color limits here
vmin = None  # e.g., vmin = 0.0
vmax = None  # e.g., vmax = 1.0

vmin = -1  # e.g., vmin = 0.0
vmax = 1  # e.g., vmax = 1.0

#################################################################

def create_video_from_frames_mp4(ts_vec, fld, frame_folder='./out_frames/', output_video='output_video.mp4', fps=10):
    """
    Create a .mp4 video from PNG frames using mp4v codec.
    """
    frame_files = [os.path.join(frame_folder, f"{fld}_frame_{n_step:08d}.png") for n_step in ts_vec]
    frame_files = [f for f in frame_files if os.path.exists(f)]
    if not frame_files:
        print("No frames found in", frame_folder)
        return

    first_frame = cv2.imread(frame_files[0])
    height, width, layers = first_frame.shape

    # Use mp4v codec for .mp4 files (widely supported)
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

    for frame_file in frame_files:
        img = cv2.imread(frame_file)
        if img is not None:
            video.write(img)
        else:
            print(f"Warning: could not read {frame_file}")

    video.release()
    print(f"Video saved as {output_video}")


def create_video_from_frames(ts_vec, fld, frame_folder='./out_frames/', output_video='output_video.mov', fps=10):
    """
    Create a .mov video from PNG frames using H.264 codec.
    """
    frame_files = [os.path.join(frame_folder, f"{fld}_frame_{n_step:08d}.png") for n_step in ts_vec]
    frame_files = [f for f in frame_files if os.path.exists(f)]
    if not frame_files:
        print("No frames found in", frame_folder)
        return

    first_frame = cv2.imread(frame_files[0])
    height, width, layers = first_frame.shape

    # Use H.264 codec for .mov files
    fourcc = cv2.VideoWriter_fourcc(*'avc1')
    video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

    for frame_file in frame_files:
        img = cv2.imread(frame_file)
        if img is not None:
            video.write(img)
        else:
            print(f"Warning: could not read {frame_file}")

    video.release()
    print(f"Video saved as {output_video}")


# Check if the output frames folder exists, create if not
if not os.path.exists('./out_frames/'):
    os.makedirs('./out_frames/')
    
# foldername = foldername + 'results/'

plt.rcParams.update({'font.size': 5})

done = 0


x = np.linspace(0, 2*np.pi, nx)
y = np.linspace(0, 2*np.pi, nx)
z = np.linspace(0, 2*np.pi, nx)

# Define the dimensions of the reshaped arrays (nvec) [y z x]
nvec = (nx, nx, nx)

# nvec = (512, 513, 512)  # Update with the actual dimensions
# nvec = (256, 257, 256)
# # nvec = (128, 129, 128)
# nvec = (0, 0, 0)

nx = nvec[2]
ny = nvec[0]
nz = nvec[1]

for n_step in ts_vec:
    file_names = []

    file_names.append(fld + '_{:08d}.dat')
    out_framename = './out_frames/' + fld + '_frame' + '_{:08d}.png'

    # Check if the output frames file exists, continue if it does not
    if not os.path.exists(out_framename.format(n_step)):

        print(f"Processing step {n_step}...")

        # Add the path string to file names
        file_names = [file_name.format(n_step) for file_name in file_names]

        # Initialize an empty list to store the data arrays
        data_arrays = []

        # Read the data from each file and reshape
        id_fnames = -1
        for file_name in file_names:
            id_fnames = id_fnames+1
            with open(foldername + file_name, 'rb') as file:
                # # Assuming the data is stored as a 1D array of floats
                # data = np.fromfile(file, dtype=np.float64)
                # data = data.reshape(nvec)*1.0
                # # data = np.flip(data, axis=1)  # Flip along the second axis (y-axis) to ensure z-axis is ascending
                # data = data.transpose((2, 1, 0)) # Permute first and third index to match the convection [x,z,y]
                # data = np.flip(data, axis = 1)

                            # Assuming the data is stored as a 1D array of floats
            # data = np.fromfile(file, dtype=np.float64)

                total_elements = np.prod(nvec)
                data = np.memmap(file, dtype=np.float64, mode='r', shape=(total_elements,))


                data = data.reshape(nvec)*1.0
                # data = np.flip(data, axis=1)  # Flip along the second axis (y-axis) to ensure z-axis is ascending
                # data = data.transpose((2, 1, 0)) # Permute first and third index to match the convection [x,z,y]
                # data = np.flip(data, axis = 1)

                dersuff = ''
                if derivative_vec[0] != 0:
                    dersuff = '_'
                    for ider in derivative_vec:
                        derivative_x, derivative_z, derivative_y = np.gradient(data, x, z, y)
                        if ider == 1:
                            data = derivative_x
                            dersuff = dersuff+'x'
                        elif ider == 2:
                            data = derivative_y
                            dersuff = dersuff+'y'
                        elif ider == 3:
                            data = derivative_z
                            dersuff = dersuff+'z'

                # define axes:
                if slice_dir == 1:
                    hor_name = 'y'
                    ver_name = 'z'
                    hor = y
                    nhor = ny
                    ver = z
                    nver = nz
                elif slice_dir == 2:
                    hor_name = 'x'
                    ver_name = 'z'
                    hor = x
                    nhor = nx
                    ver = z
                    nver = nz
                elif slice_dir == 3:
                    hor_name = 'x'
                    ver_name = 'y'
                    hor = x
                    nhor = nx
                    ver = y
                    nver = ny

                if slice_idx == -1:
                    type_name = 'average'
                    if slice_dir == 1:
                        mean_array = np.transpose(np.mean(data, axis=0))
                    elif slice_dir == 2:
                        mean_array = np.mean(data, axis=2)
                    elif slice_dir == 3:
                        mean_array = np.mean(data, axis=1)

                else:
                    type_name = 'slice'
                    if slice_dir == 1:
                        mean_array = np.transpose(data[slice_idx-1,:,:])
                    elif slice_dir == 2:
                        mean_array = data[:,:,slice_idx-1]
                    elif slice_dir == 3:
                        mean_array = data[:,slice_idx-1,:]

                data_arrays.append(mean_array)

        # Plot each array as a heat map with coordinates

        for i, array in enumerate(data_arrays):
            if array.ndim != 2:
                print(f"Invalid shape for {file_names[i]}")
                continue
            
            N = 500j
            extent = (hor.min(),hor.max(),ver.min(),ver.max())

            HOR,VER = np.meshgrid(hor, ver)
            hors,vers = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
            aaa = np.ones((len(HOR.flatten()),2))
            aaa[:,0] = HOR.flatten()
            aaa[:,1] = VER.flatten()

            plt.figure(figsize=(figsize_val, figsize_val*0.12+figsize_val))
            
            ax = plt.gca()
            # im = ax.imshow(np.flip(resampled.T, axis = 0), cmap='jet', origin='lower', extent=extent, aspect='auto')#extent=[y.min(), y.max(), z.min(), z.max()])
            # im = ax.imshow(resampled.T, cmap='jet', origin='lower', extent=extent, aspect='auto')#extent=[y.min(), y.max(), z.min(), z.max()])
            # Use pcolormesh for non-uniform grids
            X, Y = np.meshgrid(hor, ver)

            im = ax.pcolormesh(X, Y, array.T, cmap='jet', shading='gouraud', vmin=vmin, vmax=vmax)  # smoother shading

            ax.set_aspect('equal')  # <-- Add this line to set axis equal

            plt.xlabel(hor_name, fontsize=fontsize_val)
            plt.ylabel(ver_name, fontsize=fontsize_val)
            plt.xticks(fontsize=fontsize_val, rotation=0)
            plt.yticks(fontsize=fontsize_val, rotation=0)

            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("top", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax,orientation = "horizontal")
            cax.xaxis.set_ticks_position('top')
            plt.xticks(fontsize=fontsize_val, rotation=0)
            plt.title(file_names[i]+dersuff,fontsize=fontsize_val)

            # plt.show()

            # save the figure as output frame
            plt.savefig(out_framename.format(n_step), dpi=300)

            plt.close()  # <-- Add this line
 
    else:
        print(f"File {out_framename.format(n_step)} already exists. Skipping this step.")

videoname = f'./out_frames/{fld}_{ts_vec[0]}_{ts_vec[-1]}.mp4'
create_video_from_frames_mp4(ts_vec, fld, frame_folder='./out_frames/', output_video=videoname, fps=10)
# create_video_from_frames(ts_vec, fld, frame_folder='./out_frames/', output_video='./out_frames/output_video.mov', fps=10)

