
library(Seurat)
library(loomR)
library(here)
library(tidyverse)
library(ggsci)
library(cowplot)
library(viridis)
library(reshape2)
library(pheatmap)
library(colorspace)
library(scales)
library(qs)
library(cocoframer)
theme_set(theme_cowplot())

library(ComplexHeatmap)
library(future)
library(furrr)

library(data.table)

#Download from https://github.com/AllenInstitute/cocoframer
dirrds <- 'cocoframer/00.database/download_files/'
dirfun <- 'cocoframer/00.database/function_scripts/'
diraba <- 'cocoframer/00.database/aba'

run_cocoframer <- function(obj,gene.use=NULL,assay.use= 'RNA',group.by=NULL,normalize="slice"){
# Load ABA reference ------------------------------------------------------
#reference = get_ccf_annotation()
reference = readRDS(paste0(dirrds,'reference.rds'))

# Load ontology -----------------------------------------------------------
#ont = get_mba_ontology()
ont <- readRDS(paste0(dirrds,'ont.rds'))
ont_df = flatten_mba_ontology(ont)
taxon = generate_mba_taxons(ont_df)
taxon_filt = taxon
#taxon_filt = filter_mba_ontology_children(ont_df, "CTX")
taxon_filt1 = taxon_filt %>% filter(st_level>1) %>%
  mutate(acronym = make.names(acronym))

# Process ontology --------------------------------------------------------

taxon_md = list("ACA" = c("meso", "D"),
                "AI" = c("meso", "L"),
                "AOB" = c("allo", "V"),
                "AON" = c("allo", "V"),
                "AUD" = c("neo", "D"),
                "BLA" = c("amyg", "V"),
                "BMA" = c("amyg", "V"),
                "CA" = c("allo", "M"),
                "CLA" = c("meso", "L"),
                "COA" = c("amyg", "V"),
                "DG" = c("allo", "M"),
                "ECT" = c("meso", "L"),
                "ENT" = c("allo", "V"),
                "EP" = c("allo", "V"),
                "FRP" = c("neo", "D"),
                "GU" = c("neo", "D"),
                "ILA" = c("meso", "M"),
                "LA" = c("amyg", "V"),
                "MO[a-z]" = c("neo", "D"),
                "MOB" = c("allo", "V"),
                "NLOT" = c("allo", "V"),
                "ORB" = c("neo", "D"),
                "PA$" = c("amyg", "V"),
                "PAA" = c("allo", "V"),
                "PAR" = c("allo", "M"),
                "PERI" = c("meso", "L"),
                "PIR" = c("allo", "V"),
                "PL" = c("meso", "M"),
                "POST" = c("allo", "M"),
                "PRE" = c("allo", "M"),
                "ProS" = c("allo", "M"),
                "PTL" = c("neo", "D"),
                "RSP" = c("meso", "D"),
                "SS" = c("neo", "D"),
                "SUB" = c("allo", "M"),
                "TE" = c("neo", "D"),
                "TR" = c("allo", "V"),
                "TT" = c("allo", "V"),
                "VIS" = c("neo", "D"),
                "VISC" = c("meso", "L")
                )

taxon_md_df = data.frame(prefix = names(taxon_md),
                         cortical = map_chr(taxon_md, 1),
                         pallial = map_chr(taxon_md, 2),
                         stringsAsFactors = F) %>%
  mutate(prefix = paste0("^", prefix))

taxon_md_terms = map(seq_along(taxon_md_df$prefix), function(i) {
  grep(taxon_md_df$prefix[i], taxon_filt1$acronym, value=T)
})
names(taxon_md_terms) = taxon_md_df$prefix
taxon_md_terms

taxon_md_terms_u = unlist(taxon_md_terms)
taxon_terms_df = tibble(prefix = names(taxon_md_terms),
                            acronym = taxon_md_terms)  %>%
  unnest(cols=c(acronym))

taxon_md_df  = taxon_md_df %>% left_join(taxon_terms_df)

taxon_md_df = taxon_md_df %>%
  mutate(pallial = factor(pallial, levels=c("D", "L", "M", "V")),
         cortical = factor(cortical, levels=c("neo", "meso", "allo", "amyg")))

taxon_md_df = taxon_md_df %>% left_join(taxon_filt1)
taxon_md_df1 = taxon_md_df %>% filter(n_children==0)
#taxon_md_df1 <- readRDS('taxon_md_df1.rds')

# input data preparation ---------------------------------------------------------------
#source('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/huangbaoqian/ST_PRECISION_USER_huangbaoqian/cocoframer/02.github/songbird_cells/utils/go.R')
#tfs = get_tf_genes_human()
#saveRDS(tfs,'tfs.rds')
tfs = readRDS(paste0(dirrds,'tfs.rds'))

#obj_int_filt = qread('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/huangbaoqian/ST_PRECISION_USER_huangbaoqian/cocoframer/02.github/birds_seurat_qs/HVC_RA_GABA.qs')
Idents(obj) <- group.by
obj$celltype <- obj@meta.data[,group.by]
obj_int_filt <- obj
if (is.null(gene.use)) {
obj_int_filt = FindVariableFeatures(obj_int_filt, nfeatures = 6000)
var_genes = VariableFeatures(obj_int_filt)

genes_to_use = var_genes
#genes_to_use = var_genes[var_genes %in% tfs$external_gene_name]
}else{
var_genes <- gene.use
genes_to_use = var_genes
#genes_to_use = var_genes[var_genes %in% tfs$external_gene_name]
}
ave <- AverageExpression(obj_int_filt,features = genes_to_use,assays=c(assay.use))
obj_int_avg1 <- ave[[assay.use]]
obj_int_avg1 = obj_int_avg1
#obj_int_avg1 = obj_int_avg1[rownames(obj_int_avg1) %in% tfs$external_gene_name,]

# Load ABA ISH ------------------------------------------------------------
aba_data_dir <- diraba

gene_id_fname = file.path(dirrds, "gene_id.rds")
redo = F
if (redo) {
  gene_id = get_exp_gene_relationships()
  saveRDS(gene_id, gene_id_fname)
} else {
  gene_id = readRDS(gene_id_fname)
}

gene_id_filt = gene_id %>%
  mutate(gene_symbol = toupper(gene_symbol)) %>%
  filter(gene_symbol %in% genes_to_use)

ish_dir = file.path(diraba, "ish_summaries")
dir.create(ish_dir)

ishs = map(gene_id_filt$id, function(id) {
  ish_fname_cur = file.path(ish_dir, sprintf("%s.qs", id))
  print(id)
  if (! file.exists(ish_fname_cur)) {
    ish = get_aba_ish_structure_data(id)
    qsave(ish, ish_fname_cur)
  } else {
    ish = qread(ish_fname_cur)
  }
  ish
#}, .progress=T, .options=future_options(scheduling=10))
})
ishs = ishs %>%
  set_names(gene_id_filt$id) %>%
  bind_rows(.id="id")

ishs = ishs %>% left_join(gene_id_filt)

ishs1 = ishs %>%
  mutate(acronym = make.names(acronym)) %>%
  group_by(gene_symbol, atlas_id, name, acronym) %>%
  summarize(energy_mean = mean(energy))
ish_mat = acast(ishs1, gene_symbol~acronym, value.var="energy_mean", fill=0)
ish_mat = ish_mat[,colnames(ish_mat) %in% taxon_filt1$acronym]

#ish_mat = ish_mat[intersect(rownames(ish_mat), egg_l_mm$gene),]
#rownames(ish_mat) = egg_l_zf$predicted_name[match(rownames(ish_mat), egg_l_mm$gene)]

# Correlate ---------------------------------------------------------------
source(paste0(dirfun,'spatial_funcs.R'))
obj_int_avg2 = obj_int_avg1[rownames(obj_int_avg1) %in% rownames(ish_mat),]
ish_mat2 = ish_mat[rownames(ish_mat) %in% rownames(obj_int_avg2),]
mat_a = log1p(as.matrix(obj_int_avg2)) + .1
mat_b = log1p(ish_mat2[,colnames(ish_mat2) %in% taxon_md_df1$acronym]) + .1
obj_cor = region_correlate(mat_a, mat_b, method="spearman")

# Shuffle -----------------------------------------------------------------
obj_cor_shuf_fname <- 'obj_cor_shuf.rds'

redo = T
if (redo) {
  nrep = 100
  plan(multiprocess(workers=12, gc=T))
  obj_cor_shuf = future_map(1:nrep, function(i) {
    mat_a_cur = mat_a
    mat_b_cur = mat_b
    rownames(mat_a_cur) = sample(rownames(mat_a_cur))
    obj_cor = region_correlate(mat_a_cur, mat_b_cur)


    obj_cor_df = melt(obj_cor)
    colnames(obj_cor_df) =  c("celltype", "acronym", "value")
    obj_cor_df
  }) %>% bind_rows()
  qsave(obj_cor_shuf, obj_cor_shuf_fname)
} else {
  obj_cor_shuf = qread(obj_cor_shuf_fname)
}

obj_cor_shuf_stat = obj_cor_shuf %>% group_by(celltype, acronym) %>%
  summarize(value_mean = mean(value),
            value_sd = sd(value),
            value_q99 = quantile(value, .99),
            value_q95 = quantile(value, .95))

obj_cor_df = melt(obj_cor)
colnames(obj_cor_df) = c("celltype", "acronym", "value")
obj_cor_df = obj_cor_df %>% left_join(obj_cor_shuf_stat)

# Map cor to anno ---------------------------------------------------------
taxon_filt_cur = taxon_filt %>%
  mutate(acronym = make.names(acronym)) %>%
  select(acronym, id)

obj_cor_df_filt1 = obj_cor_df %>% left_join(taxon_filt_cur)
obj_cor_id = obj_cor
colnames(obj_cor_id) = taxon_filt$id[match(colnames(obj_cor), taxon_filt$acronym)]
obj_cor_id = obj_cor_id[,!is.na(colnames(obj_cor_id))]
unique_ids = as.numeric(colnames(obj_cor_id))
unique_cts = rownames(obj_cor)
#grid_anno_cors = map(seq_along(unique_cts), function(i) array(NA, dim=dim(reference))) %>%
grid_anno_cors = map(seq_along(unique_cts), function(i) array(NA, dim=dim(reference))) %>%
  set_names(unique_cts)

ref_dt = as.data.table(reference)
setkey(ref_dt, value)
for (ct in unique_cts) {
  print(ct)
  obj_cor_df_cur = obj_cor_df_filt1 %>% filter(celltype==ct)
  unique_ids = unique(obj_cor_df_cur$id)
for (i in unique_ids) {
  obj_cor_df_cur1 = obj_cor_df_cur %>% filter(id==i)
  inds = as.matrix(ref_dt[.(i), .(V1, V2, V3)])
  grid_anno_cors[[ct]][inds] = obj_cor_df_cur1$value[1]
}
}

obj_cor_df_top = obj_cor_df_filt1 %>%
  filter(value>value_q95) %>%
  group_by(celltype) %>%
  top_n(10, value)

# Calculate median positions ----------------------------------------------

ref_dt_med = ref_dt %>% group_by(value) %>%
  summarize(X = median(V1),
            Y = median(V2),
            Z = median(V3)) %>%
  rename(id = value) %>%
  left_join(taxon_filt %>% select(id, acronym))
ref_dt_med %>% filter(acronym=="RSPv1")
taxon_filt %>% filter(id==60)


# Plot ISH, dense ---------------------------------------------------------
source(paste0(dirfun,'cocoframer_functions.R'))
out_dir <- './figures/'
dir.create(out_dir)
plane = "coronal"
resolution = 10
cts = rev(unique(obj_cor_df_filt1$celltype))
slices = seq(10, 400, resolution)
save(cts,plane,reference,out_dir,slices,grid_anno_cors,taxon_filt,normalize,file='out.RData')
map(cts, function(ct) {
  print(ct)
  if (plane=="coronal") {
    height = dim(reference)[2] * .01
    width = dim(reference)[3] * .01
  }
  pdf(paste0(out_dir,'fig_',ct,'.pdf'),width,height)
  map(slices, function(slice_num) {
    print(slice_num)
  ish_plot = ish_slice_heatmap_flat(mat=grid_anno_cors[[ct]],
                                    anno = reference,
                                    taxon = taxon_filt,
                                    slice_num = slice_num,
                                    plane = plane,
                                    normalize = normalize,
                                    colorset = c(brewer_pal(palette= "Blues")(9)[c(1,9)])
  ) +
    labs(title=ct)


  ish_plot <- ish_plot + ggtitle(paste0(ct,'_',slice_num))
  print(ish_plot)

  })
  dev.off()
})
}


replot_cocoframer <- function(inrdata){
load(inrdata)
resolution = 10
slices = seq(10, 400, resolution)
map(cts, function(ct) {
  print(ct)
  if (plane=="coronal") {
    height = dim(reference)[2] * .01
    width = dim(reference)[3] * .01
  }
  pdf(paste0(out_dir,'fig_',ct,'.pdf'),width,height)
  map(slices, function(slice_num) {
    print(slice_num)
  ish_plot = ish_slice_heatmap_flat(mat=grid_anno_cors[[ct]],
                                    anno = reference,
                                    taxon = taxon_filt,
                                    slice_num = slice_num,
                                    plane = plane,
                                    normalize = normalize,
                                    colorset = c(brewer_pal(palette= "Blues")(9)[c(1,9)])
  ) +
    labs(title=ct)


  ish_plot <- ish_plot + ggtitle(paste0(ct,'_',slice_num))
  print(ish_plot)

  })
  dev.off()
})
}
