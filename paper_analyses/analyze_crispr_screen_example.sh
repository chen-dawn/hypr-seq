HYPRDIR=/seq/lincRNA/Projects/MIP/ben/201211_hypr_crispr_screen_code_update/hypr-seq/
PAPERDIR=$HYPRDIR/paper_analyses/
DATADIR=/seq/lincRNA/Projects/MIP/ben/191106_rerun_for_paper_figs/
SAMPLE=Old_EM30_Deep-1

mkdir -p output

python $PAPERDIR/analyze_crispr_screen.py --collapsed ${DATADIR}/${SAMPLE}/${SAMPLE}_collapsed_emul.count_matrix.txt \
	--uncollapsed ${DATADIR}/${SAMPLE}/${SAMPLE}_uncollapsed_emul.count_matrix.txt --barcode_stats ${DATADIR}/${SAMPLE}/tmp/${SAMPLE}.collapse_bc_stats.txt \
	--guide_info $PAPERDIR/guide_info_example.txt --target_genes GATA1 --housekeeping_genes GAPDH --outdir ./output/ --name $SAMPLE --codedir $HYPRDIR