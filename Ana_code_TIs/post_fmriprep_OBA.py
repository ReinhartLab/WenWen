#!/usr/bin/env python
# coding: utf-8

# In[1]:


# %matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
from nilearn import plotting
from nilearn.image import load_img, math_img
from nilearn.glm.first_level import FirstLevelModel
import pandas as pd
import numpy as np
import os
from os.path import join
import glob
from scipy.stats import norm
from bids import BIDSLayout


# In[2]:


# define directory
myDir = {
    'bids': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/',
    'prep': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/',
    'fs': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/freesurfer/',
    'fig': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/Figs/',
    'res': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/Results/',
    'ana': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/Ana_code/',
    'templates': '/projectnb2/viscog01/Wen/IBL_TI_fMRI/Ana_code/templates/',
}
os.environ['SUBJECTS_DIR'] = '/projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/freesurfer'

layout = BIDSLayout(myDir['bids'])
print(layout)
all_subs = layout.get_subjects()
# all_subs


# In[3]:


subj = open(join(myDir['bids'],'subj_single.txt'), 'r').readlines()
subj=subj[0].strip()

# subj = 'S02'
subj_name= [element for element in all_subs if subj in element] # full name of this subj
# print(subj_name)
subj_name = subj_name[0]
subj_dir = join(myDir['prep'], f'sub-{subj_name}',f'ses-{subj_name}')
# print(subj_dir)


# In[4]:


func_task = 'OBA'
# inspace = 'T1w' # perform glm in individual space
inspace = 'MNI152NLin2009cAsym' # perform glm in ？ space

img_name = glob.glob(join(subj_dir,'func', f'*{func_task}*{inspace}*desc-preproc_bold.nii.gz'))
# print(img_name)
fmri_img = load_img(img_name[0])
tr = fmri_img.header.get_zooms()[3]
vox_size = fmri_img.header.get_zooms()[0]
# print(tr)

from nilearn.image import mean_img
mean_img = mean_img(fmri_img)

img_name = glob.glob(join(subj_dir,'func', f'*{func_task}*{inspace}*desc-brain_mask.nii.gz'))
fmri_mask = load_img(img_name[0])


# In[5]:


# Read events
events = pd.read_csv(join(myDir['ana'],'timing',func_task,f'{subj}.csv'))
# print(events)

# motion parameters jointly observed with fMRI acquisitions

conf_file = glob.glob(join(subj_dir,'func', f'*{func_task}*desc-confounds_timeseries.tsv'))
# print(conf_file)
add_reg_names = ['white_matter','global_signal','framewise_displacement','trans_x', 'trans_y', 'trans_z','rot_x', 'rot_y', 'rot_z']  # Replace with your actual column names
# add_reg_names = ['framewise_displacement','trans_x', 'trans_y', 'trans_z','rot_x', 'rot_y', 'rot_z']  # Replace with your actual column names
confounds_glm = pd.read_csv(conf_file[0], delimiter='\t',usecols=add_reg_names)
confounds_glm = confounds_glm.fillna(0) # replace nan with zeros, mostly for the first element in FD    


# In[6]:


hrf_model = 'spm + derivative + dispersion'
# hrf_model = 'spm'
# signal_scaling: False, int or (int, int), default=0
#     If not False, fMRI signals are scaled to the mean value of scaling_axis given, which can be 0, 1 or (0, 1). 0 refers to mean scaling each voxel with respect to time, 1 refers to mean scaling each time point with respect to all voxels & (0, 1) refers to scaling with respect to voxels and time, which is known as grand mean scaling. Incompatible with standardize (standardize=False is enforced when signal_scaling is not False).

fmri_glm = FirstLevelModel(t_r=float(tr),slice_time_ref=0.5,mask_img=fmri_mask,
                          noise_model='ar1',hrf_model=hrf_model,
                          drift_model='cosine',
                          high_pass=1./100,smoothing_fwhm=2*vox_size,
                          signal_scaling=(0,1), 
                          minimize_memory=False,n_jobs=-2)
fmri_glm = fmri_glm.fit(fmri_img, events,confounds_glm)


# In[7]:


# define contrast
design_matrix = fmri_glm.design_matrices_[0]
contrasts = {
    'Image_Scramble': np.where(design_matrix.columns.str.contains('derivative|dispersion'),0,
                         np.where(design_matrix.columns.str.contains('Image'), 1,
                              np.where(design_matrix.columns.str.contains('Scrambled'), -1,0))),
    'Image_Scramble_loc1': np.where(design_matrix.columns.str.contains('derivative|dispersion'),0,
                             np.where(design_matrix.columns.str.contains('Image_loc1'), 1,
                                  np.where(design_matrix.columns.str.contains('Scrambled_loc1'), -1,0))),
    'Image_Scramble_loc2': np.where(design_matrix.columns.str.contains('derivative|dispersion'),0,
                             np.where(design_matrix.columns.str.contains('Image_loc2'), 1,
                                  np.where(design_matrix.columns.str.contains('Scrambled_loc2'), -1,0))),
    'Image_Scramble_loc3': np.where(design_matrix.columns.str.contains('derivative|dispersion'),0,
                             np.where(design_matrix.columns.str.contains('Image_loc3'), 1,
                                  np.where(design_matrix.columns.str.contains('Scrambled_loc3'), -1,0))),
}

# from nilearn.plotting import plot_contrast_matrix

# for key, values in contrasts.items():
#     plot_contrast_matrix(values, design_matrix=design_matrix)
#     plt.suptitle(key)
# plt.show()


# In[9]:


# define masks

# using OBA template from Nancy Kanwisher's GSS CVS registration
# https://web.mit.edu/bcs/nklab/GSS.shtml
from nilearn.image import index_img
img_name = join(myDir['templates'],'cvs_object_parcels', 'fROIs-fwhm_5-0.0001.nii')
oba_gss_atlas = load_img(img_name)
# print(oba_gss_atlas.shape)

# Squeeze the 4D image to 3D
oba_gss_atlas = index_img(oba_gss_atlas, 0)
# print(oba_gss_atlas.shape)


# In[ ]:


# %pip install atlasreader

from nilearn.plotting import plot_stat_map,plot_img, show
from nilearn.glm.thresholding import threshold_stats_img
from nilearn.reporting import get_clusters_table, make_glm_report
from atlasreader import create_output
from pathlib import Path

def run_glm_analysis(thresh_p, cluster_vox_num, hc,isow):
    z_maps = {}
    fig, axes = plt.subplots(len(contrasts), 1, figsize=(15, 8))

    for i, (contrast_name, contrast_values) in enumerate(contrasts.items(), 1):
        out_dir = Path(myDir['fig']) / subj / func_task 
        out_dir.mkdir(parents=True, exist_ok=True)
        
        # delete previous atlasreader output if any
        files_to_delete = glob.glob(join(out_dir,'atlasreader*'))
        for file_path in files_to_delete:
            if os.path.isfile(file_path):
                os.remove(file_path)
                
        res_filename = join(myDir['res'], f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.nii.gz')
        pic_name = join(myDir['fig'], f'{subj}_{func_task}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.png')
        if os.path.exists(res_filename) and os.path.exists(pic_name) and isow == 0:
            return
            
        tmp_zmap = fmri_glm.compute_contrast(contrast_values, output_type='z_score')
        z_maps[contrast_name] = tmp_zmap

        thresholded_map, threshold = threshold_stats_img(z_maps[contrast_name], mask_img=fmri_mask, alpha=thresh_p, height_control=hc)
        thresholded_map.to_filename(res_filename)

        # get location of significant clusters in an Atlas
        create_output(thresholded_map, cluster_extent=cluster_vox_num, voxel_thresh=norm.ppf(1 - thresh_p / 2),
                      direction='pos', outdir=out_dir)

        # rename output
        old_filepath = join(out_dir,'atlasreader_clusters.csv')
        if os.path.exists(old_filepath):
            new_filepath = join(out_dir,f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_atlasreader_clusters.csv')
            os.rename(old_filepath, new_filepath)
        old_filepath = join(out_dir,'atlasreader_peaks.csv')
        if os.path.exists(old_filepath):
            new_filepath = join(out_dir,f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_atlasreader_peaks.csv')
            os.rename(old_filepath, new_filepath)
        old_filepath = join(out_dir, 'atlasreader.png')
        if os.path.exists(old_filepath):
            new_filepath = join(myDir['fig'],f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_atlasreader.png')
            os.rename(old_filepath, new_filepath)

        out_directory = Path(myDir['fig']) / 'atlasreader'
        out_directory.mkdir(exist_ok=True)

        # Use pathlib to list and rename files
        old_dir = Path(out_dir)
        png_files = old_dir.glob('atlasreader_cluster*.png')
        if any(png_files):
            _ = [file.rename(out_directory / f"{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_{file.name}") for file in png_files]

        # export report of first-level glm
        report = make_glm_report(fmri_glm, contrasts=contrasts[contrast_name], bg_img=mean_img, threshold=threshold,
                                 alpha=thresh_p, height_control=hc)
        report.save_as_html(join(myDir['fig'], f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_GLM_report.html'))

        # # cluster
        # table = get_clusters_table(z_maps[contrast_name], stat_threshold=threshold, cluster_threshold=cluster_vox_num)
        # table.to_csv(join(myDir['fig'], f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_cluster_table.csv'))

        plot_stat_map(z_maps[contrast_name], bg_img=mean_img, threshold=threshold,
                      display_mode='z', cut_coords=3, black_bg=True,
                      title=f"{contrast_name}, fdr p<{thresh_p:.3f}, threshold={threshold:.3f}", axes=axes[i - 1])

    plt.savefig(pic_name)
    plt.show()


# In[ ]:


def plot_glm_and_mask(thresh_p, cluster_vox_num, hc,isow):  
    
    contrast_name = 'Image_Scramble' 
    res_filename = join(myDir['res'], f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}_masked.nii.gz')
    pic_name = join(myDir['fig'], f'{subj}_{func_task}_glm_and_mask_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.png')

    if os.path.exists(res_filename) and os.path.exists(pic_name) and isow == 0:
        return
        
    # Load GLM Result of all locations, set Image > Scrambled
    glm_result = load_img(join(myDir['res'], f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.nii.gz'))
    glm_result = image.math_img('np.clip(img, 0, None)', img=glm_result)
        
    # resample and binarize, mask and glm result are in MNI space, arbitary thresh of 1.96 considering value interpolate
    oba_mask_resampled = image.resample_to_img(oba_gss_atlas, glm_result)
    oba_mask_resampled = image.math_img('img > 1.96', img=oba_mask_resampled)# mask p<0.05

    # Apply mask to the GLM result
    masked_glm_result = image.math_img('glm_result * roi_mask', glm_result=glm_result, roi_mask=oba_mask_resampled)
    masked_glm_result.to_filename(res_filename)

    fig, axes = plt.subplots(3, 1, figsize=(5, 10))

    plotting.plot_roi(combined_mask_img, title='OBA Mask', axes=axes[0])
    plotting.plot_stat_map(glm_result, title='GLM Result', axes=axes[1])
    plotting.plot_stat_map(masked_glm_result, title='Masked GLM Result', cut_coords=None, display_mode='ortho', colorbar=True, axes=axes[2])

    plt.savefig(pic_name)
    plt.show()


# In[ ]:


# thresh_p_values = [0.01, 0.005, 0.0001]
# cluster_vox_num_values = [3, 5, 8, 10]
# hc_values = ['fpr', 'fdr']

thresh_p_values = [0.001]
cluster_vox_num_values = [3]
hc_values = ['fpr']
isow = 1 # 1=overwrite,0= not
for thresh_p in thresh_p_values:
    for cluster_vox_num in cluster_vox_num_values:
        for hc in hc_values:
            run_glm_analysis(thresh_p, cluster_vox_num, hc,isow)
            plot_glm_and_mask(thresh_p, cluster_vox_num, hc,isow)


# In[ ]:


from nilearn import datasets, surface
from nilearn import plotting
from nilearn.datasets import fetch_surf_fsaverage

# Load fsaverage surfaces
fsaverage = datasets.fetch_surf_fsaverage(mesh='fsaverage7')

# # Load curvature data and compute sign
curv_right = surface.load_surf_data(fsaverage.curv_right)
curv_right_sign = np.sign(curv_right)

curv_left = surface.load_surf_data(fsaverage.curv_left)
curv_left_sign = np.sign(curv_left)

viewAng = ['lateral', 'posterior']
col = len(viewAng)
for thresh_p in thresh_p_values:
    for cluster_vox_num in cluster_vox_num_values:
        for hc in hc_values:

            contrast_name = 'Image_Scramble' 
            glm_result = load_img(join(myDir['res'],f'{subj}_{func_task}_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.nii.gz'))
            glm_result = image.math_img('np.clip(img, 0, None)', img=glm_result)

            # Convert glm_result to surface texture
            texture = surface.vol_to_surf(glm_result, fsaverage.pial_right)
            textureL = surface.vol_to_surf(glm_result, fsaverage.pial_left)

            fig = plt.figure(figsize=(8, 8))

            for c, view in enumerate(viewAng, 1):
                ax = fig.add_subplot(col,2, 2*c-1, projection='3d')
                im = plotting.plot_surf_stat_map(
                    fsaverage.infl_left, textureL, hemi='left',
                    colorbar=True,vmax=5,threshold=norm.ppf(1-thresh_p/2), bg_map=curv_left, view=view, axes=ax
                )
            
                ax.set_title(f'{contrast_name} ({view.capitalize()})\nUncorrected p<{thresh_p}')
                
                ax = fig.add_subplot(col,2, c*2, projection='3d')
                plotting.plot_surf_stat_map(
                    fsaverage.infl_right, texture, hemi='right',
                    colorbar=True, vmax=5,threshold=norm.ppf(1-thresh_p/2), bg_map=curv_right, view=view, axes=ax
                )
                ax.set_title(f'{contrast_name} ({view.capitalize()})\nUncorrected p<{thresh_p}')
                
            plt.savefig(join(myDir['fig'],f'{subj}_{func_task}_surf_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.png'))
                # plt.savefig(join(myDir['fig'],f'{subj}_{func_task}_surf_{contrast_name}_{hc}_thresh{thresh_p}_clusterVoxN{cluster_vox_num}.png'),dpi=600)
            plt.show()

