source /usr/bin/activate pyscenic

inscp='pyscenic/00.scripts/'
inrds='00.rds/hs_ex.rds'
ingroup='sub'

Rscript ${inscp}/create_loom_input.R -i ${inrds} -d ${ingroup} -l out

sh ${inscp}/pyscenic_from_loom.sh -i out.loom -n 20

Rscript ${inscp}/calcRSS_by_scenic.R -l aucell.loom -m metadata_subset.xls -c ${ingroup}
