# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import numpy as np
import copy as cp
import sys, math, cmath
import time as tm
import pandas as pd
import os
import datetime
import random as rd
import matplotlib
from matplotlib import pyplot as plt
from scipy import stats
import re
import mmap
import glob
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 14


#read the data
np.random.seed(1)
import scanpy as sc
wt1 = sc.read_csv("/Users/qingbowang/Desktop/slide_seq/WT_EM_1.csv", first_column_names=True)
wt2 = sc.read_csv("/Users/qingbowang/Desktop/slide_seq/WT_EM_2.csv", first_column_names=True)
wt1.var['batch'] = "wt1"
wt2.var['batch'] = "wt2"
wt1 = wt1.T
wt2 = wt2.T
dkd1 = sc.read_csv("/Users/qingbowang/Desktop/slide_seq/DKD_EM_1.csv", first_column_names=True)
dkd2 = sc.read_csv("/Users/qingbowang/Desktop/slide_seq/DKD_EM_2.csv", first_column_names=True)
dkd1.var['batch'] = "dkd1"
dkd2.var['batch'] = "dkd2"
dkd1 = dkd1.T
dkd2 = dkd2.T

conall = wt1.concatenate(wt2, dkd1, dkd2)

conall.raw = conall #save the raw data.
conall.var.index = conall.var.index.str.lower().str.capitalize().str.replace("Atp6vo1b2","Atp6v1b2") #capital + lower, for later convenience
sc.pp.calculate_qc_metrics(conall, inplace=True,  percent_top=None) #basic QC metrics
orig_genes = conall.var.index



#check distributions per batch (supllmentary)
def wt_or_dkd(int):
    if int<2: return ("wt")
    else: return ("dkd")
conall.obs["wt_or_dkd"] = conall.obs.batch.apply(lambda x: wt_or_dkd(int(x)))

plt.hist(np.log2(conall[conall.obs["batch"]=="0"].obs.total_counts+1), bins=40, label="batch1 (wt)",
         color="#00ffff22", edgecolor='#00ffff', histtype='stepfilled')
plt.hist(np.log2(conall[conall.obs["batch"]=="1"].obs.total_counts+1), bins=40, label="batch2 (wt)",
         color="#0040ff22", edgecolor='#0040ff', histtype='stepfilled')
plt.hist(np.log2(conall[conall.obs["batch"]=="2"].obs.total_counts+1), bins=40, label="batch3 (ob)",
         color="#ff000022", edgecolor='#ff0000', histtype='stepfilled')
plt.hist(np.log2(conall[conall.obs["batch"]=="3"].obs.total_counts+1), bins=40, label="batch4 (ob)",
         color="#ffbf0022", edgecolor='#ffbf00', histtype='stepfilled')
plt.xlabel("log2(total counts + 1)")
plt.ylabel("number of cells")
mu = np.mean(np.log2(conall.obs.total_counts+1))
sigma = np.std(np.log2(conall.obs.total_counts+1))
plt.axvline(x=mu - 2*sigma, color="black")
plt.axvline(x=mu + 2*sigma, color="black")
plt.legend(loc='upper right')
matplotlib.rcParams['font.size'] = 14
plt.savefig("/Users/qingbowang/Downloads/figsx_a_n_total_cnt_cutoff.pdf", dpi=500)

#KS test
for i in range(4):
    for j in range(i+1, 4):
        print ([i,j])
        print (stats.ks_2samp(np.log2(conall[conall.obs["batch"]==str(i)].obs.total_counts+1),
                                    np.log2(conall[conall.obs["batch"]==str(j)].obs.total_counts+1)))
        print ("\n")
print (conall.shape) #checking the data size before filtering everything


#filter by 2SD
hk = 'Nm_009735_b2m'
mu_hk = np.mean(np.log2(np.array(conall[:,hk].X)+1))
sigma_hk = np.std(np.log2(np.array(conall[:,hk].X)+1))
mu = np.mean(conall.obs.n_genes_by_counts)
sigma = np.std(conall.obs.n_genes_by_counts)
mu2 = np.mean(np.log2(conall.obs.total_counts+1))
sigma2 = np.std(np.log2(conall.obs.total_counts+1))
print (conall.shape) #before filtering
sc.pp.filter_cells(conall, min_genes=np.ceil(mu-2*sigma)) #min= 6 genes
sc.pp.filter_cells(conall, max_genes=np.floor(mu+2*sigma)) #nax = 22 genes
print (conall.shape) #after num. gene filtering
sc.pp.filter_cells(conall, min_counts=np.ceil(2**(mu2-2*sigma2)-1)) #56 counts
sc.pp.filter_cells(conall, max_counts=np.floor(2**(mu2+2*sigma2)-1)) # 666 counts
print (conall.shape) #after num. counts filtering
print(sum((conall[:,"Nm_009735_b2m"].X>np.ceil(2**(mu_hk-2*sigma_hk)-1))&(conall[:,"Nm_009735_b2m"].X<2**(mu_hk+2*sigma_hk)-1))) #b2m expression filtering
conall = conall[(conall[:,"Gfp"].X==0)&(conall[:,"Gapdh_mouse"].X==0)&(conall[:,"Gapdh"].X==0)&
                (conall[:,"Nm_009735_b2m"].X>np.ceil(2**(mu_hk-2*sigma_hk)-1))&(conall[:,"Nm_009735_b2m"].X<2**(mu_hk+2*sigma_hk)-1)]
print (conall.shape) #neg.cont expression expression filtering


#normalize per cell
conall.layers['raw'] = conall.X
conall.layers['log1p'] = np.log2(conall.X+1)
sc.pp.normalize_per_cell(conall, counts_per_cell_after=1e4)
conall.layers['normed'] = conall.X #log of raw kept as layer
#and take the log (we may or may not do this for future..)
sc.pp.log1p(conall)
conall.layers['normed_log1p'] = conall.X #log of raw kept as layer
#and scale the gene variance
sc.pp.scale(conall, max_value=None, zero_center=False)
conall.layers['normed_log1p_scaled'] = conall.X #log of raw kept as layer

#reduce dimension by PCA, 20pcs
sc.tl.pca(conall, svd_solver='arpack', n_comps=20)
#sc.pl.pca_overview(conall, color='Napsa') #e.g. Napsa gene
#sc.pl.pca(conall, color='Napsa', use_raw=False) #e.g.2
#sc.pl.pca_variance_ratio(conall, log=False)

#clustering
#go with nei 100, res 0.7, the emperically best parameter
n_nei = 100
res = 0.7
sc.pp.neighbors(conall, n_neighbors=n_nei, n_pcs=6)
sc.tl.leiden(conall, resolution=res) #leiden better than louvain, the tutorial says
sc.tl.umap(conall)
#sc.pl.umap(conall, color=['leiden', 'Napsa', 'Ctgf', 'Umod', "Nphs2"], palette='tab20', legend_loc="on data",layer="normed_log1p_scaled")
#plt.savefig("/Users/qingbowang/Desktop/slide_seq/0710/umap_nei{0}_res{1}.png".format(n_nei, res)) #Just sanity check


#do plot the full result UMAP for supplementery
matplotlib.rcParams['font.size'] = 16
sc.pl.umap(conall[conall.obs.leiden!="20"], color='leiden', legend_loc="on data", title="", palette='tab20', layer="normed_log1p_scaled")
plt.savefig("/Users/qingbowang/Downloads/figsx_c_umap_full.pdf", dpi=500)#for final ver

#manually change the color to make it similar enough to the main figure:
from matplotlib import cm
tab = cm.tab20.colors #and looking at https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html manually
cls_full = [tab[2],tab[3],"tab:brown","olivedrab","tab:green",
            tab[14],"tab:red",tab[15],"tab:cyan","darkgreen",
            "lightgreen","olivedrab","tab:green","tab:olive","tab:green",
            "tab:blue","olivedrab","tab:purple","tab:brown","tab:pink"]
sc.pl.umap(conall[conall.obs.leiden!="20"], color='leiden', legend_loc="on data", title="", palette=cls_full, layer="normed_log1p_scaled")
plt.savefig("/Users/qingbowang/Downloads/figs7_c_umap_full_colormatch.pdf", dpi=500)#for final ver



#heatmap as well (will take long: ~=30 min)
conall.obs["cluster_id"] = conall.obs.leiden
sc.pl.umap(conall, palette=cls_full)
sc.pl.heatmap(conall[conall.obs.leiden!="20"], var_names = conall.var.index, groupby="cluster_id", use_raw=False, layer="normed", figsize=(8,6))
plt.savefig("/Users/qingbowang/Downloads/figs7_b_heatmap_full.pdf", dpi=400, bbox_inches = 'tight')#for final ver

#added 0216: checking the variance explained by 6pcs
#conall.uns['pca']['variance_ratio'][:6]
#conall.uns['pca']['variance_ratio'].plot.bar()


#additional QCs manually:

#1. remove GFP et al from the plot
conall_slim = conall[:,pd.Series(conall.var.index).apply(lambda x: x not in ['Gapdh', 'Gfp', 'Gapdh_mouse',"Nm_009735_b2m", "Rpl13a", "Nm_013556_hprt"])]


# remove / combine some of the clusters (manually, w/ Jamie's eye)
conall_slim = conall_slim[(conall_slim.obs.leiden!="20")]
#combine 0 and 1
conall_slim.obs.leiden = conall_slim.obs.leiden.replace("0","1")
#combine 2 and 18
conall_slim.obs.leiden = conall_slim.obs.leiden.replace("18","2")
#3, 4, and 12, and 14
#14 newly added (0802)
conall_slim.obs.leiden = conall_slim.obs.leiden.replace("4","3").replace("12","3").replace("14","3")
#11 and 16
conall_slim.obs.leiden = conall_slim.obs.leiden.replace("16","11")

#merging new cluster 2 and 4 as well, which means to merge cluster 11 and 3
conall_slim.obs.leiden = conall_slim.obs.leiden.replace("11","3")

#before removing cluster 5 and 7, make sure that the overall expression is not crazily different per clusters
sc.pl.umap(conall_slim, color=['leiden', "n_counts", "log1p_total_counts", "n_genes"], palette='tab20', legend_loc="on data",
           layer="normed_log1p_scaled")
#further remove cluster 2 and 5, and do the plots, for the intermediate visualization
conall_slim = conall_slim[(conall_slim.obs.leiden!="5")]
conall_slim = conall_slim[(conall_slim.obs.leiden!="7")]


#final number count, wt and dkd
conall.obs["wt_or_dkd"].value_counts()
conall_slim.obs["wt_or_dkd"].value_counts()


#changing the cell type ID order to let the order look nicer (make it from 0, 1, 2, ... )
celltype_mapper = pd.DataFrame(conall_slim.obs.leiden.value_counts())
celltype_mapper["new_celltype_ID"] = np.arange(celltype_mapper.shape[0])+1
celltype_mapper["tentative_celltype_ID"] = celltype_mapper["new_celltype_ID"]*100
#first, map to tentative ID, and then map to new ID
clus_tmp = conall_slim.obs.leiden.apply(lambda x: celltype_mapper.loc[x,"tentative_celltype_ID"])
clus_tmp = pd.Series(np.array(clus_tmp)) #delete this annoying ordered dictionary structure
celltype_mapper.index = np.array(celltype_mapper["tentative_celltype_ID"])
conall_slim.obs["clus_new"] = np.array(clus_tmp.apply(lambda x: celltype_mapper.loc[x,"new_celltype_ID"]))
conall_slim.obs["clus_new"] = conall_slim.obs["clus_new"].astype(str)
conall_slim.obs.clus_new.value_counts()


#ordered heatmap / dotplot etc　

###name the cells
celltype_names = {"1":"PCT", "2":"Endo", "3":"PEC", "4":"Macrophages", "5":"CD-PC", "6":"PCT Muc1+", "7":"PCT Slc22a7+", "8":"TAL/DCT", "9":"CD-IC", "10":"Mesangial", "11":"Podocyte"}
clus_new_tmp = pd.Series(np.array(conall_slim.obs.clus_new))#to keep the order... categories sucks...
conall_slim.obs["clus_names"] = np.array(clus_new_tmp.apply(lambda x: celltype_names[x]))

sc.pl.umap(conall_slim, color=['clus_names'], palette='tab20', legend_loc="on data",
           layer="normed_log1p_scaled") #sanity check
#sort the categorical order
conall_slim.obs["clus_names"] = conall_slim.obs.clus_names.cat.reorder_categories(celltype_names.values())

#sanity check again
sc.pl.umap(conall_slim, color=['clus_names'], palette='tab20', legend_loc="on data",
           layer="normed_log1p_scaled")

#another sanity check: is the expression right? (e.g. Muc1+ really Muc1+ ?) -> yes it is
sc.pl.umap(conall_slim, color=['clus_names', "Muc1", "Slc22a7"], palette='tab20', legend_loc="on data",
           layer="normed_log1p_scaled") #yes this is right


#prep for supplementary: UMAP showing wt cells in grey and DKD cells in color
celltype_and_wtdkd = pd.DataFrame(conall_slim.obs[["clus_names","wt_or_dkd"]])
celltype_and_wtdkd.clus_names = np.array(celltype_and_wtdkd.clus_names)
celltype_and_wtdkd.wt_or_dkd = np.array(celltype_and_wtdkd.wt_or_dkd)
celltype_and_wtdkd.loc[celltype_and_wtdkd.wt_or_dkd=="wt","clus_names"] = "wt"
conall_slim.obs["wt_or_dkd_cellname"] = celltype_and_wtdkd.clus_names
conall_slim.obs["wt_or_dkd_cellname"] = conall_slim.obs.wt_or_dkd_cellname.cat.reorder_categories((list(celltype_names.values()) + ["wt"])) #adding wt as another categ


#change the colors to make it look clearer compared to the original
import matplotlib
cl = matplotlib.colors.TABLEAU_COLORS
#do the tab10 + two more grean colors, specific for PCT related ones
cl #provies the dict. let's manually add some. -> get some from https://matplotlib.org/3.1.0/gallery/color/named_colors.html
cl["tab:my_dense_green"] = "darkgreen"
cl["tab:my_light_green"] = "palegreen"
conall_slim.obs.clus_names.cat.reorder_categories(["CD-IC", "CD-PC", "Endo", "Macrophages", "Mesangial", "PCT", "PCT Muc1+", "PCT Slc22a7+",  "PEC", "Podocyte", "TAL/DCT"], inplace=True)
conall_slim.obs.wt_or_dkd_cellname.cat.reorder_categories(["CD-IC", "CD-PC", "Endo", "Macrophages", "Mesangial", "PCT", "PCT Muc1+", "PCT Slc22a7+",  "PEC", "Podocyte", "TAL/DCT", "wt"], inplace=True)
palletes = [cl['tab:blue'], cl['tab:cyan'], cl['tab:orange'], cl['tab:red'], cl['tab:purple'], cl['tab:green'],
            cl['tab:my_dense_green'], cl['tab:my_light_green'], cl['tab:brown'], cl['tab:pink'], cl['tab:olive']]
matplotlib.rcParams['font.size'] = 12
#main UMAP figure:
sc.pl.umap(conall_slim, color=['clus_names'], palette=palletes,legend_loc="on data",
           layer="normed_log1p_scaled") #this is the final version
plt.savefig("/Users/qingbowang/Downloads/umap_legendonsite.pdf", bbox_inches = 'tight', dpi=500)

#wt included in different color(for supp.):
matplotlib.rcParams['font.size'] = 14
palletes = [cl['tab:blue'], cl['tab:cyan'], cl['tab:orange'], cl['tab:red'], cl['tab:purple'], cl['tab:green'],
            cl['tab:my_dense_green'], cl['tab:my_light_green'], cl['tab:brown'], cl['tab:pink'], cl['tab:olive'], '#000000ff']
sc.pl.umap(conall_slim, color=['wt_or_dkd_cellname'], palette=palletes,
           layer="normed_log1p_scaled") #this is the final version, vs wt
plt.savefig("/Users/qingbowang/Downloads/umap_vs_wt_black.pdf", bbox_inches = 'tight', dpi=500)



#dotplots:
gnms = ["Aqp1","Atp11a","Enpp2","Lrp2","Miox","Napsa","Slc34a1","Muc1", "Slc22a7", "Aqp2",  "C1qa", "C1qb", "Ctgf", "Itga8", "Ptn", "Nphs2", "Synpo", "Wt1", "Ehd3", "Kdr", "Pecam1", "Plvap", "Aqp6", "Atp6v1b2", "Umod", "Slc12a1"]
conall_slim.var.index = conall_slim.var.index.str.replace("Atp6vo1b2","Atp6v1b2")
conall_slim.obs.clus_names.cat.reorder_categories(["PCT", "PCT Muc1+", "PCT Slc22a7+",  "CD-PC", "Macrophages", "Mesangial", "PEC", "Podocyte", "Endo", "CD-IC", "TAL/DCT"], inplace=True)
gnms = ["Aqp1","Atp11a","Enpp2","Lrp2","Miox","Napsa","Slc34a1","Muc1", "Slc22a7", "Aqp2",  "C1qa", "C1qb", "Ctgf", "Itga8", "Ptn", "Wt1", "Nphs2", "Synpo", "Ehd3", "Kdr", "Pecam1", "Plvap", "Aqp6", "Atp6v1b2", "Umod", "Slc12a1"]
sc.pl.dotplot(conall_slim, var_names=gnms, groupby="clus_names", use_raw=False, layer="normed", figsize=(7,5))
plt.savefig("/Users/qingbowang/Downloads/dotplot_reordered_fontsize18.pdf", bbox_inches = 'tight', dpi = 500, transparent=True)


#heatmap:
sc.pl.heatmap(conall_slim, var_names=gnms, groupby="clus_names", use_raw=False, layer="normed", figsize=(8,6))
plt.savefig("/Users/qingbowang/Desktop/slide_seq/heatmap_full_0104.png", bbox_inches = 'tight', dpi = 500)


#corr with atlas scRNA-seq
genes = gnms
atlas = pd.read_csv("~/Desktop/slide_seq/dge_hvgs.csv",sep=",", index_col=0)
atlas = atlas.T
atlas_clus = pd.read_csv("~/Desktop/slide_seq/cell_cluster_outcome.csv",sep=",", index_col=0)
print (atlas.shape[0], atlas_clus.shape[0], len(np.intersect1d(atlas_clus.index, atlas.index))) #okay.
#subset to hypr-seq genes:
print (len(genes), len(np.intersect1d(atlas.columns, genes))) #only 17 genes...?? But fine.
atlas_for_hypr = atlas.loc[:,np.intersect1d(atlas.columns, genes)]#still only 18 genes...??
atlas_for_hypr["cluster"] = np.array(atlas_clus.cluster)
atlas_for_hypr_ave = atlas_for_hypr.groupby("cluster").agg("mean")



from scipy import stats
import seaborn as sns
hypr_raw_df = pd.DataFrame(conall_slim[:,gnms].X) #let's do not raw
hypr_raw_df.columns = gnms
hypr_raw_df.index = conall_slim.obs.clus_new
hypr_ave = hypr_raw_df.groupby("clus_new").agg("mean")
hypr_ave = hypr_ave.loc[:,atlas_for_hypr_ave.columns]
hypr_ave.sort_index(inplace=True)
hypr_ave.index = np.array(hypr_ave.index) #hate this categorical structure
corrs = pd.DataFrame(0, index=hypr_ave.index, columns=atlas_for_hypr_ave.index)
for i in corrs.index:
    for j in corrs.columns:
        corrs.loc[i,j] = stats.pearsonr(hypr_ave.loc[i,:],atlas_for_hypr_ave.loc[j,:])[0]
        print ("done {0},{1}".format(i,j))
atlas_order = [1, 2, 8, 0, 6, 4, 10, 3, 5, 15, 11, 9] #removing the hypr clusters that are not interesting: [7, 12, 13, 14, 16, 17, 18, 19]
corrs_toplot = corrs[atlas_order]
corrs_toplot.index = pd.Series(corrs_toplot.index).apply(lambda x: celltype_names[x])
atlas_celltype_names = {"1":"Glomerular Endo", "2":"Fenestrated Endo", "8":"Endo", "0":"PCT", "6":"PCT", "4":"Immune",
                        "10":"CD-PC", "3":"DCT", "5":"TAL", "15":"CD-A-IC", "11":"Mesangial", "9":"Podocyte"}
corrs_toplot.columns = pd.Series(corrs_toplot.columns.astype(str)).apply(lambda x: atlas_celltype_names[x])
corrs_toplot_correct = pd.concat([corrs_toplot.iloc[:,3], corrs_toplot.iloc[:,:3], corrs_toplot.iloc[:,5:]], axis=1) #change the order, delete one PCT, 20200616
newcolnames = list(corrs_toplot_correct.columns)
corrs_toplot_correct.columns = newcolnames
plt.figure(figsize=(6,4.5))
sns.heatmap(corrs_toplot_correct, cmap="bwr", cbar_kws={'label': 'pearson r'})
plt.xlabel("")
plt.savefig("/Users/qingbowang/Downloads/vs10x_corrected.pdf", bbox_inches = 'tight')

