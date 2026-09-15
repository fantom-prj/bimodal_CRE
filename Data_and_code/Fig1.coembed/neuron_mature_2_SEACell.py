import numpy as np
import pandas as pd
import scanpy as sc

import SEACells

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100


ad = sc.read('./scDART_n_SEACell/neuron_mature2/coembed.h5ad')

n_SEACells = round(ad.obs.apply(len)[1]/100)
build_kernel_on = 'X_umap' 
n_waypoint_eigs = 10

model = SEACells.core.SEACells(ad, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)

model.construct_kernel_matrix()
M = model.kernel_matrix
model.initialize_archetypes()


SEACells.plot.plot_initialization(ad, model, save_as="./scDART_n_SEACell/neuron_mature2/SEACells.initialization.png")

model.fit(min_iter=10, max_iter=50)

SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')
SEACell_ad.obs["old_obs_names"] = ""
for name in SEACell_ad.obs_names:
	old_names = ad.obs[ad.obs["SEACell"] == name].index.tolist()
	SEACell_ad.obs.loc[SEACell_ad.obs_names == name, 'old_obs_names'] = "\t".join(old_names)
SEACell_ad.obs.to_csv("./scDART_n_SEACell/neuron_mature2/metadata.tsv", sep = "\t", index = False)

