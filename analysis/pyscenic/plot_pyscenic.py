import gc
gc.collect()

exec(open('functions/function_pyscenic_plot.py').read())
sc.settings.njobs = 20

FIGURES_FOLDERNAME=os.getcwd()+'/figures/'
RESULTS_FOLDERNAME=os.getcwd()+'/'
Path(FIGURES_FOLDERNAME).mkdir(parents=True,exist_ok=True)
os.chdir(FIGURES_FOLDERNAME)
inauc=RESULTS_FOLDERNAME+"/aucell.loom"
inctx=RESULTS_FOLDERNAME+'/ctx.csv'

BASE_URL = "https://motifcollections.aertslab.org/v10nr_clust/logos/"

auc_mtx,lregulons,gregulons,oregulons,df_motifs,sig_regulons=get_object(inauc=inauc,inctx=inctx,BASE_URL = "https://motifcollections.aertslab.org/v10nr_clust/logos/")
df_regulons = get_logowebsite(sig_regulons,base_url = "https://motifcollections.aertslab.org/v10nr_clust/logos/")

#https://github.com/aertslab/pySCENIC/issues/360
##################################################
#Regulon activity binarization
########################################
auc_mtx.columns=gregulons
bin_mtx, thresholds = binarize(auc_mtx)
bin_mtx.to_csv('bin_mtx.csv')
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv('thresholds.csv')
#bin_mtx = pd.read_csv('bin_mtx.csv', index_col=0)
#thresholds=pd.read_csv('thresholds.csv',index_col=0).threshold
#bin_plot(intf='TEAD4(37g)',auc_mtx=auc_mtx,thresholds=thresholds)

incolum='sub'
#bin_heatmap
meta=pd.read_csv(RESULTS_FOLDERNAME+'/metadata_subset.xls',sep='\t',header=0)
meta['ingroup']=meta[incolum]
adata=sc.read(RESULTS_FOLDERNAME+'/out.loom')
adata.obs=meta
adata.obs['cell_type']=adata.obs[incolum]
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
adata.raw = adata
sc.pp.log1p(adata)
incolor= sns.color_palette("dark")+sns.color_palette("Paired")+sns.color_palette("Set2")+sns.color_palette("Set3")+sns.color_palette("pastel")+sns.color_palette("pastel6")

name_gregulons=gregulons
import copy
tmp=copy.copy(adata)
#mpl.pyplot.close()
adata=bin_heatmap(adata=adata,auc_mtx=auc_mtx,sig_regulons=sig_regulons,incolum=incolum,incolor=incolor,bin_mtx=bin_mtx,ftype='png',label='out_',name_gregulons=gregulons,hei=50,yticklabels=True)

ingroups=['ingroup','orig.ident']

df_obs=adata.obs
ntop=10
adata=plot_tsne(adata,ingroups=ingroups,incolor=incolor)
adata.write('tsne.h5ad')
df_results=plot_zscore(df_obs,ftype='png',incolum='cell_type')
rss=plot_rss_top(auc_mtx,adata,incolum='cell_type',ntop=ntop)
toprss=plot_rss_heatmap(auc_mtx,adata,incolum='cell_type',ntop=ntop,rss_cellType=rss)
#plot_1rss_dotplot(rss,incelltype='L1/2_IT',ntop=5)

top5=df_results.sort_values('Z', ascending=False).groupby(by='cell_type').head(ntop).sort_values('cell_type')
inregulons=['Regulon('+str(i)+')' for i in top5['regulon'].tolist()[0:7]]
plot_regulon_tsne(adata,inregulons,incolor,label='out')
plot_regulon_violin(auc_mtx,df_obs,sig_regulons,incolor,oregulons,gregulons,groupby='cell_type',rscore=4.0)
##################################################
#export regulons list
i=0
new_dict={}
for key, value in lregulons.items():
    new_dict[gregulons[i]]=value
    i+=1

with open('regulons_genelist.txt','w') as f1:
    for key, value in new_dict.items():
        f1.write('>'+str(key)+' : '+str(value)+'\n')

rss_t=rss.transpose()
rss_t.to_csv('regulons_rss_scores.xls',index=True, sep='\t',header=True)

df1=pd.read_csv('regulons_rss_scores.xls',sep='\t',index_col=0)
dict1={}
for incol in df1.columns:
    dict1[incol]=df1[incol].sort_values(ascending=False).head(10).index.to_list()

df1=pd.DataFrame(dict1)
df1.to_csv('regulons_rss_top10.xls',index=True, sep='\t',header=True)
