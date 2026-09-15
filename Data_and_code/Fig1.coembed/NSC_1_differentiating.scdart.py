import sys
sys.path.append('./')

import numpy as np
import pandas as pd

import torch
from sklearn.decomposition import PCA

import graph_tools as gt
import scDART.utils as utils
import scDART.TI as ti
import scDART

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

seeds = [0]
latent_dim = 8
learning_rate = 3e-4
n_epochs = 500
use_anchor = False
reg_d = 1
reg_g = 1
reg_mmd = 1
ts = [30, 50, 70]
use_potential = True

counts_rna = pd.read_csv("/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/RNA.count.csv", index_col = 0)
counts_atac = pd.read_csv("/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/ATAC.count.csv", index_col = 0)
label_rna = pd.read_csv("/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/RNA.label.txt", header = None)
label_atac = pd.read_csv("/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/ATAC.label.txt", header = None)
coarse_reg = pd.read_csv("/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/region2gene.csv", sep = ",", index_col = 0).values


scDART_op = scDART.scDART(n_epochs = n_epochs, latent_dim = latent_dim,         ts = ts, use_anchor = use_anchor, use_potential = use_potential, k = 10,         reg_d = 1, reg_g = 1, reg_mmd = 1, l_dist_type = 'kl', seed = seeds[0],         device = torch.device('cuda' if torch.cuda.is_available() else 'cpu'))

scDART_op = scDART_op.fit(rna_count = counts_rna.values, atac_count = counts_atac.values, reg = coarse_reg, rna_anchor = None, atac_anchor = None)
z_rna, z_atac = scDART_op.transform(rna_count = counts_rna.values, atac_count = counts_atac.values, rna_anchor = None, atac_anchor = None)


pca_op = PCA(n_components = 8)
z = pca_op.fit_transform(np.concatenate((z_rna, z_atac), axis = 0))
z_rna_pca = z[:z_rna.shape[0],:]
z_atac_pca = z[z_rna.shape[0]:,:]
_ = utils.plot_latent(z1 = z_rna_pca, z2 = z_atac_pca, anno1 = label_rna, 
                    anno2 = label_atac, mode = "joint", save = None, 
                    figsize = (15,7), axis_label = "PCA")
_ = utils.plot_latent(z1 = z_rna_pca, z2 = z_atac_pca, anno1 = label_rna, 
                    anno2 = label_atac, mode = "modality", save = None, 
                    figsize = (15,7), axis_label = "PCA")



root_cell = 35
# infer the trajectory backbone 
groups, mean_cluster, T = ti.backbone_inf(np.concatenate((z_rna, z_atac), axis = 0), resolution = 0.05)

        
pca_op = PCA(n_components = 2)
z = pca_op.fit_transform(np.concatenate((z_rna, z_atac), axis = 0))
z_rna_pca = z[:z_rna.shape[0],:]
z_atac_pca = z[z_rna.shape[0]:,:]    
mean_cluster = pca_op.transform(np.array(mean_cluster))
utils.plot_backbone(z_rna_pca, z_atac_pca, groups = groups, T = T, mean_cluster = mean_cluster, mode = "joint", figsize=(10,7), save = None, axis_label = "PCA")


dpt_mtx = ti.dpt(np.concatenate((z_rna, z_atac), axis = 0), n_neigh = 10)
pt_infer = dpt_mtx[root_cell, :]
pt_infer[pt_infer.argsort()] = np.arange(len(pt_infer))
pt_infer = pt_infer/np.max(pt_infer)
# for scRNA-Seq batch
pt_infer_rna = pt_infer[:z_rna.shape[0]]
# for scATAC-Seq batch
pt_infer_atac = pt_infer[z_rna.shape[0]:]
utils.plot_latent_pt(z1 = z_rna_pca, z2 = z_atac_pca, pt1 = pt_infer_rna, pt2 = pt_infer_atac, mode = "joint", save = None, figsize = (10,7), axis_label = "PCA")


df = pd.DataFrame(z_rna_pca)
df.to_csv('/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/RNA.PCA.csv')
df1 = pd.DataFrame(z_atac_pca)
df1.to_csv('/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/ATAC.PCA.csv')
df2 = pd.DataFrame(z_rna)
df2.to_csv('/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/RNA.latent.csv')
df3 = pd.DataFrame(z_atac)
df3.to_csv('/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/NSC_1_differentiating/ATAC.latent.csv')


# In[ ]:




