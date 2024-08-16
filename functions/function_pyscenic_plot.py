import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
from IPython.display import HTML, display
from pathlib import Path
import loompy as lp
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import math
from adjustText import adjust_text

incolor= sns.color_palette("dark")+sns.color_palette("Paired")+sns.color_palette("Set2")+sns.color_palette("Set3")+sns.color_palette("pastel")+sns.color_palette("pastel6")
COLORS=incolor

def get_object(inauc,inctx,BASE_URL = "https://motifcollections.aertslab.org/v10nr_clust/logos/"):
    lf = lp.connect(inauc, mode='r+', validate=False )
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    regulons = {}
    for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
        regulons[i] =  list(r[r==1].index.values)
    lf.close()
    lregulons=regulons
    my_dict=lregulons
    gregulons=[]
    for i in my_dict.keys():
        lentf=len(my_dict[i])
        nametf=re.sub(r"\(|\)|\+", "", i)
        gregulons.append(nametf+'('+str(lentf)+'g)')
    new_dict = {}
    for i, key in enumerate(my_dict.keys()):
        new_dict[gregulons[i]] = my_dict[key]
    oregulons=list(auc_mtx.columns)
    #df_motifs = load_motifs('reg.csv')
    df_motifs = load_motifs(inctx)
    #logo plot website
    sig_regulons = load_signatures(inctx)
    return auc_mtx,lregulons,gregulons,oregulons,df_motifs,sig_regulons
    
def fetch_logo(regulon, base_url = "https://motifcollections.aertslab.org/v10nr_clust/logos/"):
    for elem in regulon.context:
        if elem.endswith('.png'):
            return '{}{}'.format(base_url, elem)
    return ""

def get_logowebsite(sig_regulons,base_url = "https://motifcollections.aertslab.org/v10nr_clust/logos/"):
    df_regulons = pd.DataFrame(data=[list(map(op.attrgetter('name'), sig_regulons)),
                                     list(map(len, sig_regulons)),
                                     list(map(fetch_logo, sig_regulons,base_url))], index=['name', 'count', 'logo']).T
    df_regulons['name']=[re.sub(r"\(|\)|\+", "", x) for x in df_regulons['name']]
    df_regulons.to_csv('motif_regulon_logo_websiet_list.tsv',index=True, sep='\t',header=True)
    return df_regulons


def bin_plot(intf:str,auc_mtx,thresholds,dpi=100,label='out',ftype='png'):
    fig, (ax1) = plt.subplots(1, 1, figsize=(4, 4), dpi=dpi)
    fig.tight_layout()
    plot_binarization(auc_mtx, intf, thresholds[intf], ax=ax1)
    plt.tight_layout()
    plt.savefig(label+'_'+intf+"_bin_plot."+ftype)


def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

def bin_heatmap(adata,auc_mtx,sig_regulons,incolum,incolor,bin_mtx,name_gregulons=None,ftype='png',label='out_',
        col_cluster=False,row_cluster=True,wei=20,hei=20,yticklabels=None):
    add_scenic_metadata(adata, auc_mtx, sig_regulons)
    ##set regulons names
    if name_gregulons is None:
        name_gregulons=auc_mtx.columns.tolist()
    lindex=[i for i, name in enumerate(adata.var.columns) if 'Regulon' not in name]
    newcolums= adata.var.columns[lindex].tolist() + [ 'Regulon('+str(i)+')' for i in name_gregulons ]
    adata.var.columns=newcolums
    lindex=[i for i, name in enumerate(adata.obs.columns) if 'Regulon' not in name]
    newcolums= adata.obs.columns[lindex].tolist() + [ 'Regulon('+str(i)+')' for i in name_gregulons ]
    adata.obs.columns=newcolums

    #process
    adata.obs['cell_type']=adata.obs[incolum]
    df_metadata=adata.obs
    df_metadata=df_metadata.sort_values(by=['cell_type'],na_position='first')
    bin_mtx=bin_mtx.loc[df_metadata.index]
    incate=adata.obs.cell_type.unique().tolist()
    adata.obs.cell_type=pd.Categorical(adata.obs.cell_type,categories=incate)
    N_COLORS = len(adata.obs.cell_type.dtype.categories)
    COLORS = incolor[0:N_COLORS]
    cell_type_color_lut = dict(zip(adata.obs.cell_type.dtype.categories, COLORS))
    df_metadata['cell_id']=df_metadata.index.tolist()
    cell_id2cell_type_lut = df_metadata.set_index('cell_id').cell_type.to_dict()
    bw_palette = sns.xkcd_palette(["white", "black"])
    ##legend1
    sns.set()
    sns.set_style("whitegrid")
    palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])
    plt.savefig(label+"_legend-on_off."+ftype)
    #legend2
    sns.set()
    sns.set(font_scale=0.8)
    palplot(sns.color_palette(COLORS), adata.obs.cell_type.dtype.categories, size=1.0)
    plt.savefig(label+"_legend-cell_type_colors."+ftype)
    #heatmap
    sns.set()
    sns.set(font_scale=1.0)
    sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
    if yticklabels is None:
        yticklabels=bin_mtx.columns.tolist()
    g = sns.clustermap(bin_mtx.T, col_colors=bin_mtx.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
            cmap=bw_palette,figsize=(wei,hei),col_cluster=col_cluster,row_cluster=row_cluster,
            dendrogram_ratio=0.1,
            yticklabels=yticklabels,cbar_pos=(0.02, 0.03, 0.03, 0.02))# (left, bottom, width, height)
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_xticks([])
    g.ax_heatmap.set_xlabel('Cells')
    g.ax_heatmap.set_ylabel('Regulons')
    g.ax_col_colors.set_yticks([0.5])
    g.ax_col_colors.set_yticklabels(['Cell Type'])
    g.cax.set_visible(False)
    g.fig.savefig(os.path.join(label+"_clustermap-on_off_heatmap."+ftype), format='png')
    return adata

def plot_tsne(adata,ingroups=['sub','major_sub'],incolor=incolor):
    sc.tl.tsne(adata, use_rep='X_aucell')
    embedding_aucell_tsne = pd.DataFrame(adata.obsm['X_tsne'], columns=[['_X', '_Y']], index=adata.obs_names)
    embedding_aucell_tsne.to_csv('AUCell_tSNE_position.tsv',index=True, sep='\t',header=True)
    sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
    sc.pl.tsne(adata, color=ingroups,
               title=ingroups, ncols=3, palette=incolor,
              save='AUCell_tSNE.png')
    return adata


def plot_zscore(df_obs,ftype='png',incolum='cell_type',zscore=None):
    df_obs['cell_type']=df_obs[incolum]
    signature_column_names = list(df_obs.select_dtypes('number').columns)
    signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
    df_scores = df_obs[signature_column_names + ['cell_type']]
    df_results = ((df_scores.groupby(by='cell_type').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
    df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
    if zscore is None:
        zscore=df_results.Z.quantile(.9)
    df_results[(df_results.Z >= zscore)].sort_values('Z', ascending=False).head()
    df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= zscore].sort_values('Z', ascending=False),
                                       index='cell_type', columns='regulon', values='Z')
    #df_heatmap.drop(index='Myocyte', inplace=True) # We leave out Myocyte because many TFs are highly enriched (becuase of small number of cells).
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))
    sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray',
                        cmap="YlGnBu", annot_kws={"size": 6})
    ax1.set_ylabel('')
    fig.savefig('regulons_zscore'+str(int(zscore))+'.'+ftype,format=ftype)
    return df_results


def plot_rss_top(auc_mtx,adata,incolum='cell_type',ntop=5,rss_cellType=None):
    adata.obs['cell_type']=adata.obs[incolum]
    if rss_cellType is None:
        rss_cellType = regulon_specificity_scores(auc_mtx, adata.obs.cell_type)
    cats = adata.obs.cell_type.unique().tolist()
    ncols = math.ceil(len(cats)/5)
    fig = plt.figure(figsize=(15, ncols*4))
    for c,num in zip(cats, range(1,len(cats)+1)):
        x=rss_cellType.T[c]
        ax = fig.add_subplot(ncols,5,num)
        plot_rss(rss_cellType, c, top_n=ntop, max_n=None, ax=ax)
        ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
        for t in ax.texts:
            t.set_fontsize(12)
        ax.set_ylabel('')
        ax.set_xlabel('')
        adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
    fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
    fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
    plt.tight_layout()
    plt.rcParams.update({
        'figure.autolayout': True,
            'figure.titlesize': 'large' ,
            'axes.labelsize': 'medium',
            'axes.titlesize':'large',
            'xtick.labelsize':'medium',
            'ytick.labelsize':'medium'
            })
    plt.savefig("cellType-RSS-top"+str(ntop)+".pdf", dpi=600, bbox_inches = "tight")
    return rss_cellType


def plot_rss_heatmap(auc_mtx,adata,incolum='cell_type',ntop=5,rss_cellType=None):
    adata.obs['cell_type']=adata.obs[incolum]
    if rss_cellType is None:
        rss_cellType = regulon_specificity_scores(auc_mtx, adata.obs.cell_type)
    cats = adata.obs.cell_type.unique().tolist()
    topreg = []
    for i,c in enumerate(cats):
        topreg.extend(
            list(rss_cellType.T[c].sort_values(ascending=False)[:ntop].index)
        )
    topreg = list(set(topreg))
    auc_mtx_Z = pd.DataFrame(index=auc_mtx.index )
    for col in list(auc_mtx.columns):
        auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
    colors = sns.color_palette('bright',n_colors=len(cats))
    colorsd = dict(zip(cats, colors))
    colormap = [ colorsd[x] for x in adata.obs['cell_type'] ]
    sns.set()
    sns.set(font_scale=0.8)
    fig = palplot( colors, cats, size=1.0)
    plt.savefig("cellType-heatmap-legend-top"+str(ntop)+".pdf", dpi=600, bbox_inches = "tight")
    sns.set(font_scale=1.2)
    g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
        yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
        cmap="YlGnBu", figsize=(21,16) )
    g.cax.set_visible(True)
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.set_xlabel('')
    plt.savefig("cellType-rss_heatmap-top"+str(ntop)+".pdf", dpi=600, bbox_inches = "tight")
    return topreg

def plot_1rss_dotplot(rss,incelltype='L1/2_IT',ntop=5):
    sns.set()
    sns.set(style='whitegrid', font_scale=0.8)
    fig, ((ax1)) = plt.subplots(1, 1, figsize=(6, 6), dpi=100)
    plot_rss(rss, incelltype, ax=ax1,top_n=ntop)
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    plt.tight_layout()
    fig.savefig(('rss_'+incelltype+'.png').replace('/','-'),format='png')

def plot_regulon_tsne(adata,inregulons,incolor,label):
    nlen=len(inregulons)
    sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
    sc.pl.tsne(adata, color=['cell_type']+inregulons,
               title=['cell_type']+inregulons, ncols=4, use_raw=False,
              save='_aucell_score_tsne_regulon_'+label+'.png', palette=incolor, cmap='viridis')

def plot_regulon_violin(auc_mtx,df_obs,sig_regulons,incolor,oregulons,gregulons,groupby='cell_type',rscore=4.0):
    aucell_adata = sc.AnnData(X=auc_mtx.sort_index())
    aucell_adata.obs = df_obs
    names = list(map(op.attrgetter('name'), filter(lambda r: r.score > rscore, sig_regulons)))
    nindex=[i for i,name in enumerate(oregulons) if name in names ]
    names=list(np.array(gregulons)[nindex])
    sc.pl.stacked_violin(aucell_adata, names, groupby=groupby,color=incolor,
              save='_regulons_violin.png')


def plot_regulons_heatmap(auc_mtx,adata,incolum='cell_type',inregulons=None,c_cluster=False,prefix='out',r_cluster=False,f_type=".pdf",f_size=(18,10),tree_h=0.2,bar_pos= (0.02, 0.8, 0.05, 0.18),font_scale=1.2):
    adata.obs['cell_type']=adata.obs[incolum]
    adata1=adata[list(auc_mtx.index)].copy()
    topreg = inregulons
    nlen=len(inregulons) 
    cats = adata1.obs.cell_type.unique().tolist()
    auc_mtx_Z = pd.DataFrame(index=auc_mtx.index)
    for col in list(auc_mtx.columns):
        auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
    colors = sns.color_palette('bright',n_colors=len(cats))
    colorsd = dict(zip(cats, colors))
    colormap = [ colorsd[x] for x in adata1.obs['cell_type'] ]
    sns.set(font_scale=0.8)
    fig = palplot( colors, cats, size=1.0)
    plt.savefig(prefix+"_cellType-heatmap-legend-inregulons"+str(nlen)+".pdf", dpi=600, bbox_inches = "tight")
    sns.set(font_scale=font_scale)
    #sns.set(font_scale=font_scale,font="Arial")
    #plt.rcParams['figure.figsize'] = f_size
    #sns.set_theme(rc={'figure.figsize':(18,10)})
    g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
        yticklabels=False, xticklabels=True, vmin=-1, vmax=3, row_colors=colormap,
        cmap="YlGnBu", col_cluster=c_cluster,row_cluster=r_cluster,
        figsize=f_size,dendrogram_ratio=tree_h,cbar_pos=bar_pos)
    g.cax.set_visible(True)
    #sns.move_legend(g,"upper left",bbox_to_anchor=(.55, .45))
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.set_xlabel('')
    #plt.rcParams['figure.figsize'] = f_size
    #plt.legend(title='cell_type', loc='bottom left', labels=cats,labelcolor=colors)
    plt.savefig(prefix+"_cellType_heatmap-inregulons"+str(nlen)+f_type, dpi=600, bbox_inches = "tight")

