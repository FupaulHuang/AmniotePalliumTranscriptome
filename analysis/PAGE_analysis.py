import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './output.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures

adata = sc.read_10x_mtx('/jdfssz3/ST_STOMICS/P22Z10200N0664/1.SH_MacacaBrain/huangbaoqian/09.in_plot/pvalb_matrix')
meta=pd.read_csv('/jdfssz3/ST_STOMICS/P22Z10200N0664/1.SH_MacacaBrain/huangbaoqian/09.in_plot/metadata_20240218.xls',sep='\t',header=0)
adata.obs['second']=meta['second']

sc.pp.recipe_zheng17(adata,n_top_genes=5000)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='second', legend_loc='on data')

pl.savefig('f1.pdf')

sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='second', legend_loc='on data')

pl.savefig('f2.pdf')
sc.tl.paga(adata, groups='second')
sc.pl.paga(adata, threshold=0.03, show=False)

pl.savefig('f3.pdf')

sc.pl.paga(adata, color=['PVALB'],threshold=0.03, show=False)
pl.savefig('f4.pdf')

sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['second','PVALB','DAAM2','NR2E1'])

pl.savefig('f5.pdf')
