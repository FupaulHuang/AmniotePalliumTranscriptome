#!/bin/bash
# **********************************************************
# * Author        : HuangFubaoqian
# * Email         : huangbaoqian@genomics.cn
# * Create time   : 2022-08-31 10:31
# * Filename      : pyscenic_from_loom.sh
# * Description   : 
# **********************************************************
#default value
input_loom=out.loom
n_workers=20
#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-h|--help:\t Usage information"
exit 1
}
#get value
while getopts :i:n:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done

database='/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/huangbaoqian/ST_PRECISION_USER_huangbaoqian/SCENIC/04.database/02.human/'
tfs=${database}/hs_hgnc_tfs.txt
feather=${database}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=${database}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
pyscenic=pyscenic

# grn
$pyscenic grn \
--num_workers $n_workers \
--output grn.tsv \
--method grnboost2 \
$input_loom  $tfs

# cistarget
$pyscenic ctx \
grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
ctx.csv \
--output aucell.loom \
--num_workers $n_workers
