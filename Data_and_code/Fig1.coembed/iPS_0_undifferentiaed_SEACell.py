import numpy as np
import pandas as pd
import anndata as ad
import SEACells

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100


# --------------------------------------------------
# Read scDART latent coordinates
# rows = cells
# columns = latent dimensions
# --------------------------------------------------

latent = pd.read_csv(
    './SEACell/iPS_0_undifferentiaed/scDART_latent.csv',
    index_col=0
)

print(latent.shape)
print(latent.head())


# --------------------------------------------------
# Create minimal AnnData object
# SEACells only needs cell IDs + embedding for assignment
# --------------------------------------------------

adata = ad.AnnData(
    X=np.zeros((latent.shape[0], 1))
)

adata.obs_names = latent.index.astype(str)

# Store full scDART latent embedding
adata.obsm['X_scDART'] = latent.to_numpy()


# --------------------------------------------------
# Define number of SEACells
# target ~100 cells per metacell
# --------------------------------------------------

n_SEACells = round(adata.n_obs / 100)

print("Number of cells:", adata.n_obs)
print("Number of SEACells:", n_SEACells)


# --------------------------------------------------
# Run SEACells on full scDART latent space
# --------------------------------------------------

build_kernel_on = 'X_scDART'
n_waypoint_eigs = 10

model = SEACells.core.SEACells(
    adata,
    build_kernel_on=build_kernel_on,
    n_SEACells=n_SEACells,
    n_waypoint_eigs=n_waypoint_eigs,
    convergence_epsilon=1e-5
)

model.construct_kernel_matrix()
M = model.kernel_matrix

model.initialize_archetypes()

sc.pp.neighbors(adata, use_rep='X_scDART')
sc.tl.umap(adata)

SEACells.plot.plot_initialization(
    adata,
    model,
    plot_basis='X_umap',
    save_as='./SEACell/iPS_0_undifferentiaed/SEACells.initialization.png'
)

model.fit(
    min_iter=10,
    max_iter=50
)


# --------------------------------------------------
# Export cell-to-SEACell assignment
# --------------------------------------------------

assignment = adata.obs[['SEACell']].copy()
assignment.index.name = 'cell_barcode'

assignment.to_csv(
    './SEACell/iPS_0_undifferentiaed/SEACell_assignment.tsv',
    sep='\t'
)

