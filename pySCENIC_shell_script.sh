
##First step using pyscenic grn 
nohup pyscenic grn \
  IPLA_NK_exp.loom \
  ~/IPLA_tumor_NK_Pyscenic/ref_file/allTFs_hg38.txt \
  --output IPLA_NK_exp_adj.csv \
  --method grnboost2 \
  --sparse \
  --num_workers 8 \
  &>pyscenic.log &
  
##Second step using pyscenic ctx 
pyscenic ctx IPLA_NK_exp_adj.csv \
    ~/IPLA_tumor_NK_Pyscenic/ref_file/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname ~/IPLA_tumor_NK_Pyscenic/ref_file/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname IPLA_NK_exp.loom \
    --output IPLA_NK_merge_regulons.csv \
    --mode custom_multiprocessing \
    --num_workers 10 \
    --mask_dropouts

##Third step using pyscenic aucell to form loom document

pyscenic aucell IPLA_NK_exp.loom \
    IPLA_NK_merge_regulons.csv \
    --output IPLA_NK_merge_SCENIC.loom \
    --num_workers 3