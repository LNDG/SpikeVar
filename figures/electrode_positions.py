import nibabel as nib
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy
# load positions
figure_dir=f"/Users/kamp/PhD/spikevar/repo/SpikeVar/figures"
position_file=f"{figure_dir}/electrodePos_clean.csv"
positions = pd.read_csv(position_file, sep=";")
positions[positions=="n/a "] = np.NAN

n_patients = len(positions)
brain_regions = positions.columns[1:]

# load mni
mni_path = "/usr/local/fsl/data/standard/MNI152_T1_0.5mm.nii.gz"
mni_image = nib.load(mni_path)
aff = mni_image.affine
mni_data = mni_image.get_fdata()

# create amygdala mask
amygdala_mask = np.zeros_like(mni_data)
for i in range(n_patients):
    for region in brain_regions[:-2]: 
        if not isinstance(positions.loc[i,region], str):
            continue
        mm_coordinates = [float(c) for c in positions.loc[i, region].replace(" ","").split(",")]
        # transform coordinates in voxel space
        x,y,z = np.array([int(c) for c in nib.affines.apply_affine(np.linalg.inv(aff), mm_coordinates)])
        # add value at electrode position
        amygdala_mask[x,y,z] = 1

# create hippocampus mask
hpc_mask = np.zeros_like(mni_data)
for i in range(n_patients):
    for region in brain_regions[2:]: 
        if not isinstance(positions.loc[i,region], str):
            continue
        mm_coordinates = [float(c) for c in positions.loc[i, region].replace(" ","").split(",")]
        x,y,z = np.array([int(c) for c in nib.affines.apply_affine(np.linalg.inv(aff), mm_coordinates)])
        hpc_mask[x,y,z] = 1

# create spheres around dots
def create_sphere(size):
    x,y,z = np.mgrid[-size:size+1,-size:size+1,-size:size+1]
    mask = x**2 + y**2 + z**2 <= size**2
    return mask.astype('int')
sphere = create_sphere(4)

hpc_mask = scipy.ndimage.convolve(hpc_mask, weights = sphere)
amygdala_mask = scipy.ndimage.convolve(amygdala_mask, weights = sphere)

# binarize
hpc_mask[hpc_mask!=0]=1
amygdala_mask[amygdala_mask!=0]=1

# save as nifti
hpc_nifti = nib.Nifti1Image(hpc_mask, aff, nib.Nifti1Header())
nib.save(hpc_nifti, "hpc_mask.nii.gz")
hpc_nifti = nib.Nifti1Image(amygdala_mask, aff, nib.Nifti1Header())
nib.save(hpc_nifti, "amygdaly_mask.nii.gz")

# create transparent to red colormap
colors = [(1,0,0,c) for c in np.linspace(0,1,100)]
cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=10)
colors = [(0,0,1,c) for c in np.linspace(0,1,100)]
cmapblue = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=10)

def plot_positions(mni_slice, hpc_mask_slice, amygdala_mask_slice, ax, view=0):
    extent = (0, hpc_mask_slice.shape[0], 0, hpc_mask_slice.shape[1])
    mni_slice[np.isclose(mni_slice,0)]=np.max(mni_slice)-10
    mni_figure = ax.imshow(mni_slice.T, cmap="gray", origin="lower", extent=extent)
    mask_figure = ax.imshow(hpc_mask_slice.T, cmap=cmapred, origin="lower", extent=extent)
    mask_figure = ax.imshow(amygdala_mask_slice.T, cmap=cmapblue, origin="lower", extent=extent)
    
    # Formatting the ticks
    views = np.arange(0,3)
    x_idx = np.min(views[views!=view])
    y_idx = np.max(views[views!=view])
    xtick_transform = lambda x, pos: nib.affines.apply_affine(aff, np.repeat(x,3)).astype(int)[x_idx]
    ax.xaxis.set_major_formatter(xtick_transform)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ytick_transform = lambda x, pos: nib.affines.apply_affine(aff, np.repeat(x,3)).astype(int)[y_idx]
    ax.yaxis.set_major_formatter(ytick_transform)
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))

    # axis labels
    ax.set_ylabel("MNI coordinates (mm)")
    if view==2:
        ax.set_xlabel("MNI coordinates (mm)")

# load 1 mm brain MNI
mni_path = "/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
mni_image = nib.load(mni_path)
mni_data = mni_image.get_fdata()

fig, ax = plt.subplots(3,1, figsize=(5,10))
x_pos, y_pos, z_pos = np.array([111, 107, 55]) # coordinates in 1mm voxelspace

# saggital
mni_slice = mni_data[x_pos, :, :]
hpc_mask_slice = hpc_mask[x_pos*2, :, :]    # electrode positions are in 0.5mm voxelspace
amygdala_mask_slice = amygdala_mask[x_pos*2, :, :]
plot_positions(mni_slice, hpc_mask_slice, amygdala_mask_slice, ax[0])

# coronal
mni_slice = mni_data[:, y_pos, :]
hpc_mask_slice = hpc_mask[:, y_pos*2, :]
amygdala_mask_slice = amygdala_mask[:, y_pos*2, :]
plot_positions(mni_slice, hpc_mask_slice, amygdala_mask_slice, ax[1], view=1)

# axial
mni_slice = mni_data[:, :, z_pos]
hpc_mask_slice = hpc_mask[:, :, z_pos*2]
amygdala_mask_slice = amygdala_mask[:, :, z_pos*2]
plot_positions(mni_slice, hpc_mask_slice, amygdala_mask_slice, ax[2], view=2)

plt.savefig("electrode_positions.pdf", bbox_inches='tight')
