#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.special import comb
import os
from os import path
from glob import glob
from scipy.stats import poisson, pearsonr, spearmanr
from collections import Counter
from scipy.stats import ttest_ind, linregress
import seaborn as sns
from scipy.io import mmread
import csv
import gzip


## IMPORT MIP TOOLS FROM DIFFERENT DIR
import sys
sys.path.insert(0, '/seq/lincRNA/Projects/MIP/mip-pipeline/bin')
import mip_tools


plt.switch_backend('agg')


params = {'pdf.fonttype': '42',
          'font.sans-serif': 'Helvetica'}
matplotlib.rcParams.update(params)
plt.style.use('seaborn-ticks')

#######################
# define what samples we are going to be using

spleen_dir_hypr = sys.argv[1] # '/seq/lincRNA/Projects/MIP/ben/191106_rerun_for_paper_figs/'
spleen_dir_10x = sys.argv[2] # 'New_EM30_Deep-2'
spleen_assignments_10x = sys.argv[3] # 'New_EM30_Deep-2'

#######################
# snakemake outputs

output_plot = '/seq/lincRNA/RAP/Paper/HyPR/Figures/Supp/figSXj.pdf'
output_stats = '/seq/lincRNA/RAP/Paper/HyPR/Figures/Supp/figSXi.txt'

#######################
# defs
immune_genes = ['Cd79a','Cd19','Ighd','Fcer2a','Cd55','Mef2c','H2-Aa','Bank1','H2-Ob','Fchsd2','Zfp318','Pxk','Cr2','Plac8','Tm6sf1','Dtx1','Cybb','Ly6a','Pdia4','Dph5','Blnk','Mzb1','Apoe','Iglv1','Ighm','Jchain','Zbtb20','Igkc','Iglc3','Iglc1','Iglc2','Vpreb3','Cd79b','Spib','Ms4a1','Siglecg','Ly6d','Fam129c','Cd24a','Hsp90ab1','Cd83','Nme1','Egr3','Ncl','Gnl3','Hspa5','Il4i1','Hsp90aa1','Mif','Ifit3','Ifi27l2a','Rnf213','Ifit2','Isg15','Trim30a','Ifi214','Ifi209','Ms4a4c','Ifi203','Mki67','Hba-a1','Hist1h1b','Hbb-bt','Car2','Stmn1','Hmgb2','C1qb','C1qc','Vcam1','C1qa','Apol7c','Slc40a1','Hmox1','Ftl1','Ctsb','Cst3','Csf1r','Cd68','Marco','Retnlg','Il1b','Csf3r','Cxcr2','Mmp9','Lcn2','S100a8','S100a9','S100a6','Msrb1','Siglech','Cd209d','Runx2','Ctsl','Mpeg1','Bst2','Tcf4','Irf8','Tyrobp','Cd14','Ccr2','S100a4','Fn1','Fcer1g','Ifitm3','Ms4a6c','Lyz2','Chil3','Psap','Ms4a6b','Bcl11b','Igfbp4','Lef1','Tcf7','Trac','Trbc1','Trbc2','Ms4a4b','Dapl1','Il7r','Cd3d','Cd3g','Foxp3','Arl4c','Ccr9','Cd4','Cd8a','Cd8b1','Thy1','Ikzf2','Il2rb','Itgb1','Izumo1r','Ahnak','Ccl5','AW112010','Klrk1','Klrc2','Klre1','Gzma','Ncr1','Nkg7','Il2ra','Ccr6','Xcr1','Itgam','Itgax','Ccr7','Mcpt4','Fcer1a','Cebpa', 'Ptprc','Col1a2','C4b','Cd209a','Cd80','Cd86','H2-Ab1','H2-Eb1','Itgae']
ppif_locus = ['Anxa11', 'Plac9a', 'Ppif', 'Zcchc24', 'Zmiz1']
hk_genes = ['Abcf1','Alas1','B2m','Cltc','G6pdx','Gapdh','Gusb','Hprt','Ldha','Pgk1','Polr1b','Polr2a','Rpl13a','Rplp0','Tbp','Tubb5']
gene_order = immune_genes + ppif_locus + hk_genes

cell_order = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'Macrophage 1', 'Macrophage 2', 'Monocyte', 'pDC', 'cDC', 'CD8 T', 'Cd4 T', 'Activated CD4 T', 'NK']

cell_type_colors = dict(zip(cell_order, sns.hls_palette(len(cell_order))))

s = pd.Series(1, index=cell_order)
s.plot.pie(colors=[cell_type_colors.get(c) for c in cell_order])
plt.ylabel('')
plt.tight_layout()
plt.savefig('/seq/lincRNA/RAP/Paper/HyPR/Figures/Supp/cell_type_colors.png')
plt.close()

#######################
# load the hypr data
spleen_files = glob(path.join(spleen_dir_hypr, '*', '*_collapsed_emul.count_matrix.txt'))

counts_hypr = mip_tools.combine_count_tables(spleen_files, remove_tails=False) # quantile=0.025

metadata_hypr = pd.DataFrame(index=counts_hypr.index, columns=['probes', 'cells'])

metadata_hypr['probes'] = metadata_hypr.apply(lambda row: '2k' if '2k' in row.name else '5k', axis=1)
metadata_hypr['cells'] = metadata_hypr.apply(lambda row: 'Top2' if 'Top2' in row.name else 'All', axis=1)
metadata_hypr['downsampled'] = metadata_hypr.apply(lambda row: 'downsampled' in row.name, axis=1)

# remove probes that aren't actually included
counts_hypr = counts_hypr[counts_hypr.columns[counts_hypr.sum() > 10]]

# split into diff datasets
all_2k = counts_hypr.loc[metadata_hypr.loc[(metadata_hypr.probes=='2k')&(metadata_hypr.cells=='All')].index.tolist()]
all_5k = counts_hypr.loc[metadata_hypr.loc[(metadata_hypr.probes=='5k')&(metadata_hypr.cells=='All')].index.tolist()]
top2_2k = counts_hypr.loc[metadata_hypr.loc[(metadata_hypr.probes=='2k')&(metadata_hypr.cells=='Top2')&(metadata_hypr.downsampled==False)].index.tolist()]
top2_2k_downsampled = counts_hypr.loc[metadata_hypr.loc[(metadata_hypr.probes=='2k')&(metadata_hypr.cells=='Top2')&(metadata_hypr.downsampled==True)].index.tolist()]
top2_5k = counts_hypr.loc[metadata_hypr.loc[(metadata_hypr.probes=='5k')&(metadata_hypr.cells=='Top2')].index.tolist()]

#######################
# compare umi counts in all vs top2

# compare counts for common genes
top2_genes = ['Ighd','Fcer2a','Cr2','Plac8','Vpreb3','Cd79b','Trbc2','Lef1','Cd8b1','Cd8a','Hsp90ab1','Cd83','Ikzf2','Il2rb','C1qb','C1qc','Ifitm3','S100a4','Retnlg','Il1b','Mki67','Hba-a1','Ifit3','Ifi27l2a','Apoe','Iglv1','Gzma','Klrk1','Lef1','Trac','Siglech','Cd209d']

all_gene_avg = all_2k.mean()
top2_gene_avg = top2_2k.mean()
top2_downsampled_gene_avg = top2_2k_downsampled.mean()

# common_genes = pd.DataFrame({'all': all_gene_avg, 'top2': top2_gene_avg})
common_genes = pd.DataFrame({'all': all_gene_avg, 'top2': top2_downsampled_gene_avg})
common_genes = common_genes.loc[top2_genes+hk_genes]

slope, intercept, r_value, p_value, std_err = linregress(common_genes['top2'], common_genes['all'])

fig, ax = plt.subplots(figsize=(5,5))

common_genes.plot.scatter('top2', 'all', ax=ax, c='k')

ax.text(1, 28, 'R = {:.2f}'.format(r_value), fontsize=12)
ax.text(1, 20, 'y = {:.2f}x + {:.2f}'.format(slope, intercept), fontsize=12)
ax.set_xlim([-1, 36])
ax.set_ylim([-1, 36])
ax.set_aspect('equal', adjustable='datalim')
ax.set_xticks([0, 5, 10, 15, 20, 25, 30, 35])
ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35])

# plt.tight_layout()
plt.savefig(output_plot)
plt.close()

umis_per_cell_all = all_2k.sum(axis=1).mean()
umis_per_cell_top2 = top2_2k.sum(axis=1).mean()
umis_per_cell_top2downsampled = top2_2k_downsampled.sum(axis=1).mean()

# decide on sample to proceed with
spleen_hypr = all_2k.copy()
spleen_hypr_assignments = metadata_hypr.loc[(metadata_hypr.probes=='2k')&(metadata_hypr.cells=='All')].copy()

#######################
# get 10x data

# load 10x data
def read_10x(pathin):
    """Return Pandas Dataframe containing 10x dataset """
    
    mat=mmread(os.path.join(pathin, "matrix.mtx.gz"))
    genes_path = path.join(pathin, "features.tsv.gz")
    gene_ids = [row[0] for row in csv.reader(gzip.open(genes_path, mode='rt'), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(gzip.open(genes_path, mode='rt'), delimiter="\t")]

    gene_final = gene_names # [x+'_'+y for x,y in zip(gene_ids,gene_names)]
    
    barcodes_path = path.join(pathin, "barcodes.tsv.gz")
    barcodes = [row[0][0:16] for row in csv.reader(gzip.open(barcodes_path, mode='rt'), delimiter="\t")]
    
    DGE=pd.DataFrame(mat.toarray())
    
    DGE.index=gene_final
    DGE.columns=barcodes
    
    return DGE.T

spleen_10x = read_10x(spleen_dir_10x)
assignments_10x = pd.read_table(spleen_assignments_10x, index_col=0)


#######################
# assign hypr

# build pseudobulk profiles
def create_pseudo_bulk(counts_df, clusters_df, cluster_col, genes=gene_order):
    pseudo_bulk = {}
    for cluster, df in clusters_df.groupby(cluster_col):
        bcs = df.index.tolist()
        pseudo_bulk[cluster] = counts_df.loc[bcs, genes].mean()
    return pd.DataFrame(pseudo_bulk) 

rescaled_normalized_10x = mip_tools.column_normalize_matrix(mip_tools.rescale_matrix(spleen_10x)[gene_order])
pseudo_bulks_10x_normalized = create_pseudo_bulk(rescaled_normalized_10x, assignments_10x, 'cluster')

# assign hypr cells
def assign_cells_by_corr(counts_df, cell_bc, pseudo_bulks, genes=gene_order):
    pseudo_bulks_copy = pseudo_bulks.copy()
    pseudo_bulks_copy['cell'] = counts_df.loc[cell_bc, genes]
    corrs = pd.DataFrame(pseudo_bulks_copy).corr()['cell'].drop('cell')
    return corrs.idxmax()

rescaled_normalized_hypr = mip_tools.column_normalize_matrix(mip_tools.rescale_matrix(spleen_hypr))

spleen_hypr_assignments['cluster'] = spleen_hypr_assignments.apply(lambda row: assign_cells_by_corr(rescaled_normalized_hypr, row.name, pseudo_bulks_10x_normalized), axis=1)
pseudo_bulks_hypr_normalized = create_pseudo_bulk(rescaled_normalized_hypr, spleen_hypr_assignments, 'cluster')


#####
# gene expression heatmaps

# plot heatmaps of 10x and hypr data
def plot_heatmap(counts_df, assignments_df, assignment_col, cell_order=cell_order, gene_order=gene_order, cell_type_colors=cell_type_colors, col_cluster=False, outfig=None):
    bc_order = []
    
    for c in cell_order:
        bc_order += assignments_df.loc[assignments_df[assignment_col]==c].index.tolist()
        
    counts_reordered = counts_df.loc[bc_order, gene_order]
    
    row_colors = assignments_df.loc[bc_order, assignment_col].map(cell_type_colors)
    
    g = sns.clustermap(counts_reordered, row_cluster=False, col_cluster=col_cluster, figsize=(10,10), row_colors=row_colors)

    # do some fiddling with the axes
    g.ax_row_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_xlim([0,0])
    
    if outfig:
        g.savefig(outfig)
    
    plt.close()



plot_heatmap(rescaled_normalized_10x, assignments_10x, assignment_col='cluster', outfig=output_plot.replace('.pdf', '.10x_clustermap.pdf'))
plot_heatmap(rescaled_normalized_hypr, spleen_hypr_assignments, assignment_col='correlation_cluster', outfig=output_plot.replace('.pdf', '.hypr_clustermap.pdf'))

top2_gene_new_order = ['Cd19','Cd79a','Ighd','Fcer2a','Cr2','Vpreb3','Mki67','Nme1','Ighm','C1qb','C1qc','Ifitm3','Itgam','Retnlg','Il1b','Siglech','Itgax','Trac','Cd8a','Cd4','Foxp3','Ikzf2','Gzma','Ncr1'] 

plot_heatmap(rescaled_normalized_10x, assignments_10x, assignment_col='cluster', gene_order=top2_gene_new_order, col_cluster=False, outfig=output_plot.replace('.pdf', '.10x_clustermap.top2.pdf'))
plot_heatmap(rescaled_normalized_hypr, spleen_hypr_assignments, assignment_col='cluster', gene_order=top2_gene_new_order, col_cluster=False, outfig=output_plot.replace('.pdf', '.hypr_clustermap.top2.pdf'))


#####
# cell type frequencies barplots

# compare cell-type frequencies
freqs_hypr = spleen_hypr_assignments['cluster'].value_counts(normalize=True)
freqs_10x = assignments_10x['cluster'].value_counts(normalize=True)

freqs = pd.DataFrame({'hypr': freqs_hypr, '10x': freqs_10x})

freqs_melted = freqs.rename_axis('Cell Type').reset_index().melt(id_vars='Cell Type', value_vars=['10x', 'hypr'], var_name='Assay', value_name='Fraction of cells')
order = freqs.sort_values('10x', ascending=False).index.tolist()

fig, ax = plt.subplots()

sns.barplot(x='Cell Type', y='Fraction of cells', hue='Assay', data=freqs_melted, order=order, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

plt.tight_layout()
plt.savefig(output_plot.replace('.pdf', '.freqs.pdf'))
plt.close()

cell_corr = pearsonr(freqs['hypr'], freqs['10x'])[0]

#####
# correlation matrix of cell type expression profiles

def correlate_two_dfs(df1, df2, cell_order):
    return df1.apply(lambda s: df2.corrwith(s)).loc[cell_order,cell_order]

def plot_corr_heatmap(df1, df2, df1_name, df2_name, cell_order=cell_order, savefig=None):
    corrs = correlate_two_dfs(df1, df2, cell_order)
    fig, ax = plt.subplots()
    sns.heatmap(corrs, vmin=0, vmax=1, square=True, ax=ax)
    ax.set_xlabel(df1_name)
    ax.set_ylabel(df2_name)
    plt.tight_layout()
    if savefig:
        plt.savefig(savefig)
    return corrs

corrs = plot_corr_heatmap(pseudo_bulks_10x_normalized, pseudo_bulks_hypr_normalized, '10x', 'HyPR', savefig=output_plot.replace('.pdf', '.celltype_corr.pdf'))

gene_corr = np.mean(np.diagonal(corrs))

#####
# scatter all genes
gene_expression_hypr = spleen_hypr[gene_order].mean()
gene_expression_10x = spleen_10x[gene_order].mean()

expression = pd.DataFrame({'hypr': gene_expression_hypr, '10x': gene_expression_10x})
expression['log_hypr'] = np.log10(expression['hypr'])
expression['log_10x'] = np.log10(expression['10x'])

fig, ax = plt.subplots()

expression.plot.scatter('log_10x', 'log_hypr', ax=ax)

print(expression)

slope_, intercept_, r_value_, p_value_, std_err_ = linregress(expression.replace([-np.inf, np.inf], np.nan).dropna()['log_10x'], expression.replace([-np.inf, np.inf], np.nan).dropna()['log_hypr'])

plt.tight_layout()
plt.savefig(output_plot.replace('.pdf', '.gene_ex_scatter.pdf'))
plt.close()

spleen_files_probes = glob(path.join(spleen_dir_hypr, '*', '*_collapsed_emul.count_matrix_probes.txt'))

counts_hypr_probes = mip_tools.combine_count_tables(spleen_files_probes, remove_tails=False) # quantile=0.025

num_probes = sum([len(counts_hypr_probes.filter(like=g).columns) for g in top2_genes+hk_genes])


# make the probes figure:
counts_hypr_probes = counts_hypr_probes.loc[spleen_hypr_assignments.index.tolist()]

def plot_increase_in_probes(genes, cell_types, counts, assignments, threshold=0, savefig=None):
    assert len(genes) == len(cell_types)
    
    res = []
    
    for idx in range(len(genes)):
        g, cts = genes[idx], cell_types[idx]
        cells = assignments.cluster.isin(cts)
        n_cells = sum(cells)
        rest_cells = len(cells) - n_cells
        all_probes = counts.loc[cells].filter(like=g).sum().sort_values(ascending=False).index.tolist()
        best_probe = all_probes[0]
        best_probe_exp = counts.loc[cells,best_probe].mean()
        best_probe_frac = (counts.loc[cells,best_probe]>threshold).sum() / n_cells
        all_probes_exp = counts.loc[cells,all_probes].sum(axis=1).mean()
        all_probes_frac = (counts.loc[cells,all_probes].sum(axis=1)>threshold).sum() / n_cells
        rest_exp = counts.loc[~cells,all_probes].sum(axis=1).mean()
        rest_frac = (counts.loc[~cells,all_probes].sum(axis=1)>threshold).sum() / rest_cells
        
        res.append({'gene': g, 
                    'cell_types': ','.join(cts), 
                    'best': best_probe_frac, 
                    'best_exp': best_probe_exp, 
                    'all': all_probes_frac, 
                    'all_exp': all_probes_exp, 
                    'rest': rest_frac,
                    'rest_exp': rest_exp})
        
    df = pd.DataFrame(res)
    
    df['fc'] = df['all'] / df['best']
    
    df1 = df.melt(id_vars=['gene','cell_types'], value_name='frac', var_name='how')

    fig, ax = plt.subplots(figsize=(5,5))

    sns.barplot(x='gene', y='frac', hue='how', data=df1, hue_order=['best','all','rest'], ax=ax)

    ax.set_ylim(0,1)

    if savefig:
        plt.tight_layout()
        plt.savefig(savefig)
    plt.close()
    
    return df
    
df = plot_increase_in_probes(['Cd209a','Ms4a4b','Klrc2'], [['cDC','pDC'],['Cd4 T','CD8 T','Activated CD4 T'],['NK']], counts_hypr_probes, spleen_hypr_assignments, savefig=output_plot.replace('.pdf', '.probes.pdf'))


with open(output_stats, 'w') as out_stats:
    out_stats.write('line equation\ty={:.3f}x+{:.3f}\n'.format(slope, intercept))
    out_stats.write('R value\t{}\n'.format(r_value))
    out_stats.write('Num probes\t{}\n'.format(1017))
    out_stats.write('Num genes\t{}\n'.format(179))
    out_stats.write('line equation\ty={:.3f}x+{:.3f}\n'.format(slope_, intercept_))
    out_stats.write('R value\t{}\n'.format(r_value_))
    out_stats.write('HyPR cells\t{}\n'.format(len(spleen_hypr)))
    out_stats.write('10x cells\t{}\n'.format(len(spleen_10x)))
    out_stats.write('Common genes\t{}\n'.format(len(common_genes)))
    out_stats.write('Common probes\t{}\n'.format(num_probes))
    out_stats.write('All probes cells\t{}\n'.format(len(all_2k)))
    out_stats.write('Top2 probes cells\t{}\n'.format(len(top2_2k)))
    out_stats.write('Cell corr\t{}\n'.format(cell_corr))
    out_stats.write('Gene corr\t{}\n'.format(gene_corr))
    out_stats.write('Best fc\t{}\n'.format(df['fc'].max()))
    out_stats.write('Worst fc\t{}\n'.format(df['fc'].min()))
    out_stats.write('UMIs/cell all\t{}\n'.format(umis_per_cell_all))
    out_stats.write('UMIs/cell top2\t{}\n'.format(umis_per_cell_top2))
    out_stats.write('UMIs/cell top2 downsampled\t{}\n'.format(umis_per_cell_top2downsampled))
