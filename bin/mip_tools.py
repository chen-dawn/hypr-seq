import itertools
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import distance
plt.switch_backend('agg')
from scipy.stats import zscore, entropy, norm, sem, ttest_ind
from scipy.cluster import hierarchy as hc
from scipy import spatial as sp
from os import path
from glob import glob
from scipy.optimize import minimize_scalar
from sklearn.mixture import GaussianMixture

from outlier_aware_hist import *

## FILL THIS IN
acc_number_to_gene = {}
mip_barcodes = ['brainbar1', 'brainbar2', 'brainbar3', 
    'brainbar4', 'brainbar5', 'brainbar6', 'brainbar7', 'brainbar8', 
    'barcode01', 'barcode02', 'barcode04', 'barcode05', 'barcode06', 
    'barcode09', 'barcode10', 'barcode11', 'barcode12', 
    'barcode13', 'barcode14', 'barcode15', 'barcode17', 
    'barcode18', 'barcode20', 'barcode22', 'barcode23', 
    'barcode24', 'barcode25', 'barcode26', 'barcode27', 
    'barcode28', 'barcode30', 'Brainbar1', 'Brainbar2', 
    'Brainbar3', 'Brainbar4', 'Brainbar5', 
    'Brainbar6', 'Brainbar7', 'Brainbar8', 
    'SplintRbar1', 'SplintRbar2', 'SplintRbar4', 
    'SplintRbar5', 'SplintRbar6', 'SplintRbar9', 
    'SplintRbar10', 'SplintRbar11', 'SplintRbar12', 
    'SplintRbar13', 'SplintRbar14', 'SplintRbar15', 
    'SplintRbar17', 'SplintRbar18', 'SplintRbar20', 
    'SplintRbar22', 'SplintRbar23', 'SplintRbar24', 
    'SplintRbar25', 'SplintRbar26', 'SplintRbar27', 'SplintRbar28', 'SplintRbar30', 
    'Brainbar001', 'Brainbar002', 'Brainbar003', 'Brainbar004', 'Brainbar005',
    'Brainbar006', 'Brainbar007', 'Brainbar008', 'SplintRbar001', 'SplintRbar002']

nbins = 100
bc_thresh = 0.7 # minimum value for barcode assignment (# highest bc count / # total barcode counts)
# chastity = 6
min_barcode_umis = 10

# get cell barcode from long read name like NB501148:265:HYMKFBGXY:3:22608:25684:19729_TTTTCCTAATGA_CTCTCGCAGC
def extract_cbc(long_name):
    return long_name.split('_')[1]

# get base probe name from long name (remove 3P_P2 part)
def extract_base_probe_name(long_name):
    if 'prime' in long_name:
        # probe is of format NM_0012897453primeprobe1
        return long_name.split('prime')[0][:-1]
    else:
        # probe is of format NM_0124233P_P1 or HDAC6Intron4_P1
        basename = long_name.rsplit('_', 1)[0]
        if basename[-1] in ['P', 'p'] and basename[-2] != 'F':
            # probe is of form NM_0124233P_P1
            return basename[:-2]
        else:
            # probe is of form HDAC6Intron4_P1 or GFP_11/12
            return basename

# add a (umi, count) pair to the list of all such pairs for a given cell
# cells_to_umis always defined as {cbc: [(umi, count), ...], ...}
def add_read(cbc, umi, count, cells_to_umis):
    cell = cells_to_umis.get(cbc, [])

    if cell:
        cell.append((umi, count))
    else:
        cells_to_umis[cbc] = [(umi, count)]

# determines whether a given cell barcode is bad
def is_pathological(cbc):
    gc = (cbc.count('G') + cbc.count('C')) / len(cbc)
    containsN = cbc.count('N') > 0
    longestRun = max(sum(1 for _ in l) for n, l in itertools.groupby(cbc))

    if gc > 0.75 or gc < 0.25:
        return True
    elif containsN:
        return True
    elif longestRun > 5:
        return True
    elif 'ATCTCGTATG' in cbc: # checks for no CBC (or truncated one) that is actually just P7 end
        return True 
    else:
        return False

# removes all bad cell barcodes from the collection
def drop_pathological(cells_to_umis):
    return { cell: umi_coll for (cell, umi_coll) in cells_to_umis.items() if not is_pathological(cell) }

# drops cells from the dictionary that have fewer than read_threshold reads
def drop_low_cells(cells_to_umis, umi_threshold):
    return { cell: umi_coll for (cell, umi_coll) in cells_to_umis.items() if len(umi_coll) > umi_threshold }
    # return { cell: umi_coll for (cell, umi_coll) in cells_to_umis.items() if sum([x[1] for x in umi_coll]) > read_threshold and cell.count('G') < 10 }
    # remove the all/mostly G barcodes, maybe this will help??

# calculates similarity between umi sets for two cell barcodes (number between 0 and 1)
# computes ratio of cardinality of intersection of umi sets by cardinality of smaller set
def umi_similarity(cbc1, cbc2, cells_to_umis):
    umis1 = {x[0] for x in cells_to_umis[cbc1]}
    umis2 = {x[0] for x in cells_to_umis[cbc2]}

    return len(umis1 & umis2) / min(len(umis1), len(umis2))

# returns cardinality of overlap of umis between two cells
def overlap_umis(cbc1, cbc2, cells_to_umis):
    umis1 = {x[0] for x in cells_to_umis[cbc1]}
    umis2 = {x[0] for x in cells_to_umis[cbc2]}

    return len(umis1 & umis2)

# drops cells from the dictionary that have high overlap in their umis (>1 dropseq bead)
# can probably do something to speed this up, maybe don't do the testing if either of the cbcs is already in cbcs_to_remove??
def filter_multiple_beads(cells_to_umis, similarity_threshold):
    cbcs_to_remove = []

    for (cbc1, cbc2) in itertools.combinations(cells_to_umis, 2):
        similarity = umi_similarity(cbc1, cbc2, cells_to_umis)
        if similarity > similarity_threshold:
            cbcs_to_remove.append(cbc1)
            cbcs_to_remove.append(cbc2)

    return { cell: umi_coll for (cell, umi_coll) in cells_to_umis.items() if cell not in cbcs_to_remove }

# assign mip barcodes to cell barcodes
# for each row, we extract only the barcode reads and take the argmax of the columns (mip_bcs) if the fraction of counts for that bc passes a threshold 
def assign_barcode(count_matrix, chastity=6.):
    barcode_counts = count_matrix[[bc for bc in mip_barcodes if bc in count_matrix.columns]]
    argmax_barcodes = barcode_counts.idxmax(axis=1)

    for cell in count_matrix.index:
        bc1, bc2 = barcode_counts.loc[cell].nlargest(2)

        bc2 = max(1.0, bc2) # make sure bc2 non-zero
        ratio = bc1 / bc2

        if ratio > chastity:
            count_matrix.loc[cell, 'Barcode'] = argmax_barcodes.loc[cell]
        else:
            count_matrix.loc[cell, 'Barcode'] = 'Unassigned'

        # if bc1 < min_barcode_umis:
        #     # not enough barcode umis
        #     count_matrix.loc[cell, 'Barcode'] = 'Unassigned'
        # else:
        #     # enough barcode umis
        #     bc2 = max(1.0, bc2) # make sure that bc2 is not 0
        #     ratio = bc1 / bc2

        #     if ratio > chastity:
        #         count_matrix.loc[cell, 'Barcode'] = argmax_barcodes.loc[cell]
        #     else:
        #         count_matrix.loc[cell, 'Barcode'] = 'Doublet'

    # return count_matrix.Barcode.isin(['Unassigned']).sum() / len(count_matrix), count_matrix.Barcode.isin(['Doublet']).sum() / len(count_matrix)
    return count_matrix.Barcode.isin(['Unassigned']).sum() / len(count_matrix)

# pass in a dictionary of real -> error cbcs from the first whitelist, 
# now we can augment what we're writing with the errors from the first step
def get_neighboring_barcodes(old_whitelist):
    # TODO
    return None

### DEPRECATED
def write_whitelist_old(whitelist_file, in_barcode_groups, n_clusters, cbcs, labels, cells_to_umis):
    cbc_clusters = [] 

    for i in range(n_clusters):
        cluster = []

        for j in range(len(cbcs)):
            if labels[j] == i:
                cluster.append(cbcs[j])

        cbc_clusters.append(cluster)

    #print(cbc_clusters[:2])

    cells_to_umis_final = {}

    for c in cbc_clusters:
        error_bcs = []
        for bc in c:
            error_bcs += in_barcode_groups[bc]
        prototype = c.pop()
        whitelist_file.write('{}\t'.format(prototype))
        # c += error_bcs
        if c + error_bcs:
            whitelist_file.write('{}\t'.format(','.join(c + error_bcs)))
        whitelist_file.write('{}'.format(len(cells_to_umis[prototype])))
        if c + error_bcs:
            whitelist_file.write('\t{}'.format(','.join(map(str, [len(cells_to_umis[b]) for b in c] + [1 for _ in error_bcs]))))
        whitelist_file.write('\n')

        # turn this into dictionary so we can recreate final thing to plot
        umi_counts = {umi: count for (umi, count) in cells_to_umis[prototype]}

        for b in c:
            umi_coll = cells_to_umis[b]
            for (umi, count) in umi_coll:
                umi_counts[umi] = umi_counts.get(umi, 0) + count

        cells_to_umis_final[prototype] = list(umi_counts.items())

    return cells_to_umis_final

# returns the number of reads associated with the specified cell dict
def count_reads(cells_to_umis):
    return sum(list(itertools.chain(*[[x[1] for x in umi_coll] for umi_coll in cells_to_umis.values()])))

def histogram_wrapper(data, output_plot, xlabel, ylabel, title, log=True):
    data = np.asarray(data)
    outlier_aware_hist(data, *calculate_bounds(data), cdf=False)
    if log:
        plt.yscale('log', nonpositive='clip')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

def cdf_wrapper(data, output_plot, xlabel, ylabel, title):
    data = np.asarray(data)
    outlier_aware_hist(np.asarray(data), *calculate_bounds(data), cdf=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

# plots histogram and cdf of reads / umi for the given dict in output_plot and returns the average
def plot_reads_per_umi_bulk(reads_per_umi, output_plot):
    histogram_wrapper(reads_per_umi, output_plot, 'Reads / UMI', 'Number of UMIs', 'Reads / UMI')
    cdf_wrapper(reads_per_umi, output_plot, 'Reads / UMI', 'Fraction of UMIs', 'Reads / UMI')

# plots histogram and cdf of reads / umi for the given dict in output_plot and returns the average
def plot_umis_per_cell_analyze_count(umis_per_cell, output_plot):
    histogram_wrapper(umis_per_cell, output_plot, 'UMIs / cell', 'Number of cells', 'UMIs / cell')
    cdf_wrapper(umis_per_cell, output_plot, 'UMIs / cell', 'Fraction of cells', 'UMIs / cell')

# plots histogram and cdf of reads / umi for the given dict in output_plot and returns the average
def plot_histogram_and_cdf(data, output_plot, xlabel=None, ylabel=None, title=None):
    histogram_wrapper(data, output_plot, xlabel, ylabel, title)
    cdf_wrapper(data, output_plot, xlabel, ylabel, title)

# plots histogram and cdf of reads / umi for the given dict in output_plot and returns the average
def plot_umi_distance(edit_dist_per_umi, output_plot):
    plt.hist(edit_dist_per_umi, bins=11)
    plt.yscale('log', nonpositive='clip')
    plt.xlabel('Edit Distance between Sequenced UMI and Representative UMI')
    plt.ylabel('Number of UMIs')
    plt.title('Edit Distance')
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

    plt.hist(edit_dist_per_umi, cumulative=True, density=True, histtype='step', bins=11)
    plt.xlabel('Edit Distance between Sequenced UMI and Representative UMI')
    plt.ylabel('Fraction of UMIs')
    plt.title('Edit Distance')
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

# plots histogram and cdf of reads / umi for the given dict in output_plot and returns the average
def plot_reads_per_umi(cells_to_umis, output_plot, title):
    reads_per_umi = np.asarray(list(itertools.chain(*[[x[1] for x in umi_coll] for umi_coll in cells_to_umis.values()]))) # is this terrible?

    histogram_wrapper(reads_per_umi, output_plot, 'Reads / UMI', 'Number of UMIs', '{}: Reads / UMI'.format(title))
    cdf_wrapper(reads_per_umi, output_plot, 'Reads / UMI', 'Fraction of UMIs', '{}: Reads / UMI'.format(title))

    return reads_per_umi.mean()

# plots histogram and cdf of reads / cell for the given dict in output_plot and returns the average
def plot_reads_per_cell(cells_to_umis, output_plot, title):
    reads_per_cell = np.asarray([sum([x[1] for x in umi_coll]) for umi_coll in cells_to_umis.values()])

    histogram_wrapper(reads_per_cell, output_plot, 'Reads / Cell', 'Number of Cells', '{}: Reads / Cell'.format(title))
    cdf_wrapper(reads_per_cell, output_plot, 'Reads / Cell', 'Fraction of Cells', '{}: Reads / Cell'.format(title))

    return reads_per_cell.mean()

def plot_umis_per_cell(cells_to_umis, output_plot, title):
    umis_per_cell = np.asarray([len(umi_coll) for umi_coll in cells_to_umis.values()])

    histogram_wrapper(umis_per_cell, output_plot, 'UMIs / Cell', 'Number of Cells', '{}: UMIs / Cell'.format(title))
    cdf_wrapper(umis_per_cell, output_plot, 'UMIs / Cell', 'Fraction of Cells', '{}: UMIs / Cell'.format(title))

    return umis_per_cell.mean()

# plot histogram (and cdf?) of probes per cell for a given probe
def plot_probe(count_matrix, probe, output_plot, title='Raw'):
    probe_counts = count_matrix[probe]
    probe_name = acc_number_to_gene.get(probe, probe)

    histogram_wrapper(probe_counts, output_plot, 'MIPs / Cell', 'Number of Cells', '{} / Cell ({})'.format(probe_name, title))
    cdf_wrapper(probe_counts, output_plot, 'MIPs / Cell', 'Fraction of Cells', '{} / Cell ({})'.format(probe_name, title))

# plot scatterplot of probes / cell for two different probes
def plot_probe_pair(count_matrix, probe_pair, output_plot):
    probe1 = probe_pair[0]
    probe2 = probe_pair[1]

    probe1_counts = count_matrix[probe1]
    probe2_counts = count_matrix[probe2]

    probe1_name = acc_number_to_gene.get(probe1, probe1)
    probe2_name = acc_number_to_gene.get(probe2, probe2)

    plt.scatter(probe1_counts, probe2_counts)
    plt.xlabel('{}'.format(probe1_name))
    plt.ylabel('{}'.format(probe2_name))
    plt.title('{} vs {}'.format(probe1_name, probe2_name))
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

def plot_counts_per_gene(collapsed_counts, plots, gene_counts, excel_file):
    counts_per_gene = collapsed_counts.sum()

    counts_per_gene = counts_per_gene.loc[counts_per_gene > 1]

    counts_per_gene.to_csv(gene_counts, sep='\t')
    pd.DataFrame(counts_per_gene).reset_index().to_excel(excel_file, header=False, index=False)

    counts_per_gene.plot.bar()
    plt.yscale('log', nonpositive='clip')
    plt.ylabel('Counts across all cells')
    plt.xlabel('Gene')
    plt.title('Sum of genes across cells')
    plt.tight_layout()
    plots.savefig()
    plt.close()

def grab_mapping_freq_from_bowtie(log_file):
    with open(log_file) as log:
        return log.readlines()[1].strip().split()[-1].strip('(').strip(')')

def extract_probe_number(long_name):
    return long_name.split('_')[-1]

def filter_cells(unfiltered_counts):
    # unfiltered_counts is cells x probes, at this point
    # filtered = unfiltered_counts

    # collect pathological barcodes
    pathological_cbcs = [cbc for cbc in unfiltered_counts.index if is_pathological(cbc)]

    # drop cells with too many or two few UMIs (probably determined by quantiles?)

    # return a subset of the rows
    return unfiltered_counts.loc[~unfiltered_counts.index.isin(pathological_cbcs)]


    #   # drop pathological cell barcodes
#   cells_to_umis_drop_pathological = mip_tools.drop_pathological(cells_to_umis)

#   # drop all cells with too few UMIs
#   cells_to_umis_drop_low = mip_tools.drop_low_cells(cells_to_umis_drop_pathological, umi_threshold)

# gets the denominators to normalize the overlap for (used below)
def get_denom(cbc1, cbc2, sets, method):
    if method == 'geometric':
        geom_mean = np.sqrt(len(sets[cbc1]) * len(sets[cbc2]))
        return geom_mean, geom_mean
    elif method == 'raw':
        return 1, 1
    elif method == 'fraction':
        return len(sets[cbc1]), len(sets[cbc2])
    else:
        return None, None

# build and return an adjacency matrix from a groups file, not thresholded yet
# can pass in choice of method:
#   geometric: each pair is fraction of overlap / geometric mean of total umis
#   raw: each pair gets raw overlap umis
#   fraction: each CBC gets fraction of overlap to total, not symmetric 
def build_adjacency_matrix(groups, method='geometric'):
    cbcs = sorted(groups.CBC.unique())
    n_beads = len(cbcs)

    # initialize empty adjacency matrix for umi overlap
    adjacency_matrix = np.zeros([n_beads, n_beads])

    sets = {}

    # create dictionary of CBC->{UMIs}
    for c, b in groups.groupby('CBC')['UMI_Transcript']:
        sets[c] = set(b.values)

    # fill in adjacency matrix 
    for i in range(n_beads):
        for j in range(i+1,n_beads):
            cbc1 = cbcs[i]
            cbc2 = cbcs[j]
            overlap = len(sets[cbc1] & sets[cbc2])
            denom_i, denom_j = get_denom(cbc1, cbc2, sets, method=method)
            adjacency_matrix[i][j] = overlap / denom_i
            adjacency_matrix[j][i] = overlap / denom_j

    return adjacency_matrix, np.array(cbcs), n_beads, sets

# helper function for plotting heatmaps
def heatmap_helper(adjacency_matrix, labels, plots, log, clustermap_png, title='Raw'):
    # plot the heatmap
    plt.imshow(adjacency_matrix)
    plt.title('Overlap matrix: {}'.format(title))
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # plot the clustermap
    # note that this sometimes fails when there isn't a valid clustering, so we will wrap in a try block
    try:
        adjacency_matrix = pd.DataFrame(adjacency_matrix) # need labels for doing clustering stuff, so gotta get it as df

        # use seaborn to get a color mapping for each unique label for the clusters
        labels_series = pd.Series(labels)
        lut = dict(zip(labels_series.unique(), sns.hls_palette(len(set(labels_series)))))
        row_colors = labels_series.map(lut)

        # get the linkage for clustering (convert to affinity matrix, then compute linkage)
        np.fill_diagonal(adjacency_matrix.values, 1) # adjust for distance computations
        affinity_matrix = 1 - adjacency_matrix
        linkage = hc.linkage(sp.distance.squareform(affinity_matrix), method='average')
        np.fill_diagonal(adjacency_matrix.values, 0) # have to adjust back for plotting

        # plot the adjacency matrix
        cg = sns.clustermap(adjacency_matrix, row_colors=row_colors, col_colors=row_colors,
                                              row_linkage=linkage, col_linkage=linkage)

        # do some fiddling with the axes
        cg.ax_row_dendrogram.set_visible(False)
        cg.ax_row_dendrogram.set_xlim([0,0])
        cg.ax_col_dendrogram.set_visible(False)
        # cg.ax_col_dendrogram.set_xlim([0,0])
        plt.title('Overlap matrix: {} (clustered)'.format(title))
        plt.tight_layout()
        plt.savefig(clustermap_png)
        # plots.savefig()
    except Exception as e:
        log.write("Couldn't cluster adjacency matrix: {}\n".format(title))
        log.write(str(e)+'\n')
    finally:
        plt.close()


# plot adjacency values per pair of cbcs
### REMOVE VALUES ALONG DIAGONAL?
def plot_adjacency_matrix(adjacency_matrix, labels, plots, thresh, log, histogram_file, clustered_file):
    # plot the histogram
    plt.hist(np.clip(adjacency_matrix.ravel(), a_min=0, a_max=1), bins=100) # bins='auto'
    plt.axvline(thresh, label='Chosen threshold', c='r', ls='--')
    plt.yscale('log', nonpositive='clip')
    plt.title('Overlap values per pair of CBCs')
    plt.ylabel('Number of barcode pairs')
    plt.xlabel('Overlap metric')
    plt.legend()
    plt.tight_layout()
    plt.savefig(histogram_file)
    plots.savefig()
    plt.close()

    # plot the heatmap
    heatmap_helper(adjacency_matrix, labels, plots, log, clustered_file)

    # plt.imshow(adjacency_matrix)
    # plt.title('Overlap matrix')
    # plt.tight_layout()
    # plots.savefig()
    # plt.close()

    # # plot the clustermap
    # # note that this sometimes fails when there isn't a valid clustering, so we will wrap in a try block
    # try:
    #     cg = sns.clustermap(adjacency_matrix)
    #     cg.ax_row_dendrogram.set_visible(False)
    #     cg.ax_row_dendrogram.set_xlim([0,0])
    #     cg.ax_col_dendrogram.set_visible(False)
    #     # cg.ax_col_dendrogram.set_xlim([0,0])
    #     plt.title('Overlap matrix (clustered)')
    #     plt.tight_layout()
    #     plots.savefig()
    # except:
    #     log.write("Couldn't cluster adjacency matrix")
    # finally:
    #     plt.close()

    # plot the binarized matrix and cluster map (if possible)
    adjacency_matrix_binarized = binarize_adjacency_matrix(adjacency_matrix, thresh)
    heatmap_helper(adjacency_matrix_binarized, labels, plots, log, clustered_file.replace('raw', 'binarized'), title='Binarized')

# convert real valued adjacency matrix into binary yes/no connection
def binarize_adjacency_matrix(adjacency_matrix, thresh):
    return adjacency_matrix > thresh

# run connected components on the binarized adjacency matrix, return the cluster labels and the number of clusters
def do_connected_components(adjacency_matrix_filtered):
    cbc_overlap = csr_matrix(adjacency_matrix_filtered)
    connected_comp = connected_components(cbc_overlap, directed=False)
    n_clusters, labels = connected_comp[0], connected_comp[1]
    return n_clusters, labels

def plot_barcodes_per_droplet(labels, plots):
    bcs_per_drop = np.unique(labels, return_counts=True)[1]
    plt.bar(*np.unique(bcs_per_drop.clip(max=6), return_counts=True))
    plt.xlabel('BCs per droplet')
    plt.ylabel('Number of droplets')
    plt.xticks(range(1,7), ['1','2','3','4','5','>5'])
    plt.tight_layout()
    plots.savefig()
    plt.close()

def plot_barcode_stats(bc_stats, plots):
    # plot pie chart of BC percentages
    bc_stats[['Pathological', 'Singleton', 'Multiplet', 'Giant']].sum().plot.pie()
    plt.title('Fraction of barcode classes')
    plt.ylabel("")
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # plot beads per droplet
    bcs_per_drop = bc_stats.drop_duplicates('ID')['n_beads'].values
    plt.bar(*np.unique(bcs_per_drop.clip(max=6), return_counts=True))
    plt.xlabel('Beads per droplet')
    plt.ylabel('Number of droplets')
    plt.xticks(range(1,7), ['1','2','3','4','5','>5'])
    plt.tight_layout()
    plots.savefig()
    plt.close()

# objective function to minimize to get lambda from zero-trucated poisson
# x_bar is the sample mean sum over x of (fraction of drops with x beads * x)
# lambda_hat is the thing we vary to estimate the normal poisson's lambda
def zero_trunc_poisson_objective(lambda_hat, x_bar):
    if lambda_hat == 0:
        return 0
    return np.abs(x_bar - lambda_hat / (1 - np.exp(-1 * lambda_hat)))

# uses observed counts of beads per droplet to estimate the poisson rate for beads
# note -- we can't do normal MLE because we don't observe how many drops with 0 barcodes
def estimate_lambda(bc_stats):
    # bc_stats_to_include = bc_stats.loc[(bc_stats.Pathological == False) & (bc_stats.Giant == False)]
    bc_stats_to_include = bc_stats.loc[(bc_stats.Pathological == False) & (bc_stats.n_beads <= 5)] # going to ignore things with >5 beads (giants?)

    # do unique, return_counts
    # have to drop duplicates (on barcode ID) so we don't overcount 
    n_beads, n_drops = np.unique(bc_stats_to_include.drop_duplicates('ID')['n_beads'].values, return_counts=True)

    # compute observed sample mean
    x_bar = np.sum(n_beads * n_drops / np.sum(n_drops))

    optim_res = minimize_scalar(lambda x: zero_trunc_poisson_objective(x, x_bar))

    if optim_res.success:
        return optim_res.x
    else:
        return -1

def write_whitelist(whitelist, whitelist_file, sets):
    with open(whitelist_file, 'w') as wl:
        for cbcs in whitelist:
            cbc = cbcs.pop(0)
            count = len(sets[cbc])
            counts = None
            if cbcs:
                counts = [len(sets[c]) for c in cbcs]
                wl.write('\t'.join([cbc, ','.join(cbcs), str(count), ','.join(map(str, counts))])+'\n')
            else:
                wl.write('\t'.join([cbc, '', str(count), ''])+'\n')

def get_n_overlapping_umis_from_df(groups_df, cbc1, cbc2):
    groups_df['UMI_Transcript'] = groups_df['final_umi'] + '_' + groups_df['contig']
    umis_1 = groups_df.loc[groups_df.CBC == cbc1, 'UMI_Transcript']


# takes in a groups table and outputs estimates of how many barcodes per drop based on shared umis
def cluster_beads_by_droplet(groups, umi_overlap_threshold):
    # grab list of cell barcodes
    cbcs = sorted(groups.CBC.unique().values)
    n_beads = len(cbcs)

    # initialize empty adjacency matrix for umi overlap
    adjacency_matrix = np.zeros([n_beads, n_beads])

    # fill in adjacency matrix 
    for i in range(n_beads):
        for j in range(i+1,n_beads):
            cbc1 = cbcs[i]
            cbc2 = cbcs[j]
            overlap = get_n_overlapping_umis_from_df(groups, cbc1, cbc2)
            # overlap = mip_tools.overlap_umis(cbc1, cbc2, cells_to_umis_drop_low)
            adjacency_matrix[i][j] = overlap
            adjacency_matrix[j][i] = overlap

    # threshold gives binary matrix of connections (are two CBCs from the same droplet or not)
    adjacency_matrix_filt = adjacency_matrix > umi_overlap_threshold

    # turn this matrix into a sparse matrix and run connected components analysis on it
    cbc_overlap = csr_matrix(adjacency_matrix_filt)
    connected_comp = connected_components(cbc_overlap, directed=False)

    # extract the number of clusters and their labels (for cbcs in each droplet)
    n_clusters, labels = connected_comp[0], connected_comp[1]

    # # read input whitelist and generate groups of cbcs:
    # in_barcode_groups = {}

    # for cbc_line in in_whitelist:
    #   cbc_split = cbc_line.split()
    #   correct_bc = cbc_split[0]
    #   if len(cbc_split) == 4:
    #       error_bcs = cbc_split[1].split(',')
    #   else:
    #       error_bcs = []
    #   in_barcode_groups[correct_bc] = error_bcs

    # # write out list of acceptable cell barcodes and the things that should be collapsed to them
    # cells_to_umis_final = mip_tools.write_whitelist(out_whitelist, in_barcode_groups, n_clusters, cbcs, labels, cells_to_umis_drop_low)

# writes the counts per probe to excel file for jamie
def write_counts_per_probe_excel(counts, outfile):
    counts_per_probe = pd.DataFrame(counts.sum(), columns=['count']).reset_index()
    counts_per_probe['probe_num'] = counts_per_probe.apply(lambda row: int(row.gene.split('_')[-1]), axis=1)
    counts_per_probe['gene_name'] = counts_per_probe.apply(lambda row: extract_base_probe_name(row.gene), axis=1)
    to_write = counts_per_probe.groupby(['gene_name', 'probe_num']).agg(np.sum).reset_index().pivot(index='gene_name', columns='probe_num', values='count').fillna(0)
    # to_write = pd.DataFrame(counts.sum(), columns=['count']).reset_index()
    to_write.reset_index().to_excel(outfile, header=False, index=False)

def calculate_zscores(counts, outfile, excel_file, log):
    zscores = zscore(counts, axis=0)
    zscores = pd.DataFrame(zscores, columns=counts.columns, index=counts.index)
    zscores.to_csv(outfile, sep='\t', header=True, index=True)
    if excel_file:
        try:
            zscores.T.reset_index().to_excel(excel_file, header=False, index=False)
        except:
            log.write("Too many CBCs, can't write excel file in this direction, try transposed\n")
            zscores.to_excel(excel_file)

def compute_entropy(counts_file):
    return counts_file[[c for c in mip_barcodes if c in counts_file.columns]].apply(lambda row: entropy(row, base=2), axis=1)

# convert HCR count table into table for comparison with sc-geno
def read_hcr_table(path_to_table):
    ct = pd.read_table(path_to_table, index_col=0)
    return ct[[c for c in ct.columns if c in mip_barcodes]]

# read geno sam file into df for comparison with hypr
def read_geno_sam(path_to_sam):
    geno_names = ['Read', 'Flags', 'Barcode']
    geno = pd.read_table(path_to_sam, header=None, usecols=list(range(3)), names=geno_names)
    return geno

# pivot sam df into count table 
def geno_sam_to_count_table(geno_sam):
    geno_sam['cbc'] = geno_sam.apply(lambda row: row.Read.split('_')[-2],axis=1)
    geno = pd.DataFrame(geno_sam.groupby('cbc')['Barcode'].value_counts())
    geno.index = geno.index.rename(names=['cbc', 'bc'])
    ct = geno.reset_index().pivot(index='cbc', columns='bc', values='Barcode').fillna(0)
    return ct.drop('*', axis=1)

# saves a knee plot from a groups file (needs columns for CBC and UMI_Transcript)
def make_knee_plot_from_groups(groups_file, plots, cbc='CBC', umi_transcript='UMI_Transcript', method='umis'):
    if method == 'umis':
        agg_func = lambda x: len(pd.unique(x))
    elif method == 'reads':
        agg_func = len
    else:
        agg_func = None

    plt.plot(sorted(groups_file.groupby(cbc)[umi_transcript].agg(agg_func).values, reverse=True))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('CBC Rank')
    plt.ylabel(method)
    plt.tight_layout()
    plots.savefig()
    plt.close()

def argparse_boolean(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

def argparse_bc_stats(s):
    if s == 'None':
        return None
    else:
        return s

# list of probes known to be problematic
probe_exclusion_list = ['GATA1Intron5P_31', 'GAT1Int5P_5']

def get_probes_for_gene(counts, gene, probe_exclusion_list=probe_exclusion_list):
    return counts[[col for col in counts.columns if extract_base_probe_name(col) == gene and col not in probe_exclusion_list]]

plot_scale = 3

def plot_probes_per_transcript(counts, probe, plots):
    # counts is cells x probes (all probes, not just collapsed ones)
    # get all columns that match the regex
    if probe in mip_barcodes:
        subset = counts.filter(regex='^{}5P_[0-9]+'.format(probe), axis=1)
    else:
        subset = get_probes_for_gene(counts, probe)
        # subset = counts.filter(like=probe, axis=1)

    subset_sum = pd.DataFrame(subset.sum(axis=0), columns=[probe])
    subset_sum['ProbeNum'] = subset_sum.apply(lambda row: extract_probe_number(row.name), axis=1)

    subset_sum.plot.bar(x='ProbeNum', y=probe, legend=None)

    plt.title('Counts per probe: {}'.format(probe))
    plt.yscale('log', nonpositive='clip')
    plt.ylim(1, plot_scale * subset_sum[probe].values.max())
    plt.ylabel('Sum of reads per probe across cells')
    plt.tight_layout()
    plots.savefig()
    plt.close()

def plot_pairplot_of_individual_probes(counts, probe, plots, title='Raw'):
    # counts is cells x probes (all probes, not just collapsed ones)
    # get all columns that match the regex
    if probe in mip_barcodes:
        subset = counts.filter(regex='^{}5P_[0-9]+'.format(probe), axis=1)
    else:
        subset = get_probes_for_gene(counts, probe)
        # subset = counts.filter(like=probe, axis=1)

    pp = sns.pairplot(subset)

    pp.fig.suptitle('Pairplot per probe: {} ({})'.format(probe, title), y=1.01)
    plt.tight_layout()
    plots.savefig()
    plt.close()

# fit a GMM to the logged counts data for each cell for the crop-hypr experiments
def fit_and_plot_model(bc_counts, plots, stats, log, n_components=2):
    # compute logged counts
    all_logs = []

    for i in range(len(bc_counts)):
        bc_cts = bc_counts.iloc[i].sort_values(ascending=False)
        best = bc_cts[0]
        bc_cts = bc_cts / best
        for l in bc_cts.loc[(bc_cts > 0) & (bc_cts < 1)]: 
            # maybe we want to include the 1s too?
            all_logs.append(np.log(l))
            
    all_logs = np.array(all_logs)

    # fit GMM to data
    model = GaussianMixture(n_components).fit(all_logs.reshape(-1,1))

    # plot results 
    plt.hist(all_logs, bins='auto')
    plt.title('Logged count ratios')
    plt.xlabel('log(ct_i/ct_0)')
    plt.ylabel('Number of observed barcodes')
    plt.tight_layout()
    plots.savefig()
    plt.close()

    x_range = np.linspace(np.min(all_logs), 0, 1000)
    pdf = np.exp(model.score_samples(x_range.reshape(-1, 1)))
    responsibilities = model.predict_proba(x_range.reshape(-1, 1))
    pdf_individual = responsibilities * pdf[:, np.newaxis]

    plt.hist(all_logs, bins='auto', density=True, histtype='stepfilled', alpha=0.5)
    plt.plot(x_range, pdf, '-k', label='Mixture')
    plt.plot(x_range, pdf_individual, '--k', label='Components')
    plt.title('Gaussian Mixture Model Fit')
    plt.xlabel('log(ct_i/ct_0)')
    plt.ylabel('Density')
    plt.legend()
    plt.tight_layout()
    plots.savefig()
    plt.close()

    stats.write('n_components\t{}\n'.format(n_components))

    # return model
    return model

# prior probability for observing t barcodes (note that t is an index and is exclusive, so t=3 means lst[:3] = [0,1,2])
def barcode_number_prior(t, max_t=7):
    t = int(min(t, max_t))

    flat_prior = [1 for _ in range(max_t + 1)]
    flat_prior[0] = 0 # can't have 0 barcodes 

    three_bcs_prior = [0 for _ in range(max_t + 1)]
    three_bcs_prior[3] = 1 # only mass on 3 

    normal_prior = [0 for _ in range(max_t + 1)]
    normal_prior[1] = 5
    normal_prior[2] = 10
    normal_prior[3] = 20
    normal_prior[4] = 10
    normal_prior[5] = 5
    normal_prior[6] = 10

    prior = np.array(normal_prior)
    prior = prior / np.sum(prior)

    return prior[t]

# returns an unnormalized probability for each threshold (to differentiate real from background barcodes)
# takes in a prior and a model for assigning probabilities to classes for a given log-count
def compute_probability(log_transformed_counts, threshold, model, real_idx):
    p_prior = barcode_number_prior(threshold)
    
    post_probs = model.predict_proba(log_transformed_counts.values.reshape(-1, 1))
    
    p_real = 1
    p_background = 1
    
    for i in range(len(log_transformed_counts)):
        if i < threshold:
            # real
            p_real *= post_probs[i, real_idx]
        else:
            # background
            p_background *= (1 - post_probs[i, real_idx])
            
    return p_prior * p_real * p_background

def construct_bc_string(top_n_bcs):
    return '-'.join(sorted(map(standardize_probe_names, top_n_bcs)))

def extract_real_crop_barcodes(counts, model, real_cluster, plots, stats, log):
    real_barcode_assignments = pd.DataFrame(index=counts.index, columns=['Threshold', 'First Choice', 'Plus 1', 'Minus 1'])

    # compute the likelihood function for all of the barcodes 
    for bc in counts.index:
        bc_cts = counts.loc[bc].sort_values(ascending=False)
        best = bc_cts[0]
        bc_cts = bc_cts / best
        bc_cts = bc_cts.loc[bc_cts>0]
        bc_cts = np.log(bc_cts)
        if len(bc_cts):
            # for each threshold,compute the probability of seeing the data, then choose the optimal one
            t = 1 + np.argmax([compute_probability(bc_cts, x, model, real_cluster) for x in range(1, 1+int(len(bc_cts)))])
            bc_string = construct_bc_string(bc_cts.index.tolist()[:t])
            bc_string_plus_one = construct_bc_string(bc_cts.index[:t+1])
            bc_string_minus_one = construct_bc_string(bc_cts.index[:t-1])
            real_barcode_assignments.loc[bc] = [t, bc_string, bc_string_plus_one, bc_string_minus_one]
        else:
            t = 0
            real_barcode_assignments.loc[bc] = [t, '', '', '']
    
    # plot a sampling of the cells
    n_rows, n_cols = (5, 4)
    fig, ax = plt.subplots(n_rows, n_cols) #, figsize=(10,15))
    axs = ax.ravel()

    # have to add in a check here for whether we have <20 cells
    to_replace = len(counts) < n_rows * n_cols

    for idx, bc in enumerate(counts.sample(n_rows * n_cols, replace=to_replace).index):
        n_bcs = 10
        bc_cts = counts.loc[bc].sort_values(ascending=False)
        best = bc_cts[0]
        bc_cts = bc_cts / best
        axs[idx].bar(range(n_bcs), bc_cts[:n_bcs])
        axs[idx].axvline(real_barcode_assignments.loc[bc, 'Threshold'])
        axs[idx].axis('off')

    fig.suptitle('Selection of barcode thresholds', y=1)
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # plot a barplot of barcodes per droplet
    real_barcode_assignments.Threshold.value_counts().reset_index().sort_values('index').plot.bar(x='index', y='Threshold', color='blue', legend=None)
    plt.xlabel('Observed Barcodes per Droplet')
    plt.ylabel('Number of Droplets')
    plt.title('Barcodes per Droplet')
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # estimate the number of "cells" per droplet
    real_barcode_assignments['Cells'] = (real_barcode_assignments['Threshold'] + 2) // 3

    # plot a barplot of cells per droplet
    real_barcode_assignments.Cells.value_counts().reset_index().sort_values('index').plot.bar(x='index', y='Cells', color='blue', legend=None)
    plt.xlabel('Estimated Cells per Droplet')
    plt.ylabel('Number of Droplets')
    plt.title('Cells per Droplet')
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # use the distribution of cells per droplet to estimate the lambda 
    vcs = real_barcode_assignments['Cells'].value_counts().reset_index()
    x_bar = (vcs['index'] * vcs['Cells'] / vcs['Cells'].sum()).sum()
    optim_res = minimize_scalar(lambda x: zero_trunc_poisson_objective(x, x_bar))

    if optim_res.success:
        cell_lambda = optim_res.x
    else:
        cell_lambda = -1

    stats.write('avg_bc_per_droplet\t{}\n'.format(real_barcode_assignments['Threshold'].mean()))
    stats.write('estimated_cell_lambda\t{}\n'.format(cell_lambda))
    stats.write('putative_singletons\t{}\n'.format((real_barcode_assignments['Threshold']<4).mean()))
    stats.write('putative_doublets\t{}\n'.format((real_barcode_assignments['Threshold']>3).mean()))

    return real_barcode_assignments
    
# used as a helper for function below, takes in a row from the observed barcodes df (with many columns, from the 
# procedure above) and a whole dictionary for bc->guide (from plasmid) and assigns a guide (called with .apply())
def get_matching_guide_from_two_dfs(observed_row, dict_as_df, log):
    observed_bcs = observed_row['First Choice']

    match = dict_as_df.loc[dict_as_df['sorted_bc_string'] == observed_bcs, 'guide']

    if match.empty:
        return 'No guide found'
    else:
        # if len(match.values) > 1:
        #     log.write(match)
        return match.values[0]

# this helper is called a couple of times in the function below
# basically, it just takes a bc_string from a set of observed barcodes and determines whether
# the plasmid corresponding to it is good or not
# note that we are only calling this when we have found only one match in the dictionary
def check_unambiguous_bc_combo(obs_str, matches):
    assert len(matches) == 1

    matches = matches.iloc[0] # just get a series from the row
    
    barcode_passes_plasmid_qc = False
    barcode_with_single_guide = False
    guide = ''
    
    # check whether the plasmid passed QC (based on having enough reads and being in the array)
    if matches['sufficient_reads'] and matches['real']:
        barcode_passes_plasmid_qc = True
        # check whether this is a unique guide or if there are multiple consistent guides
        if matches['single_guide']:
            barcode_with_single_guide = True
            guide = matches['guide']
    
    return barcode_passes_plasmid_qc, barcode_with_single_guide, guide
    
# called in a loop below, this function takes a row from the "observed barcodes" table (for one crop-hypr cbc)
# and attempts to assign it to a guide
# it does a lot of if statements which probably could be simplified, but seems to be doing something reasonable
# for now
# returns a dictionary of flags to be used to compute stats, and also a guide (if there was one, otherwise 'None')
def attempt_to_assign_one_barcode(row, dictionary):
    observed_str = row['First Choice']
    
    observed_in_dict = False
    observed_barcodes_unique_plasmid = False
    had_to_filter_low_abundance = False
    ambiguous_barcode_combination = False
    barcode_passes_plasmid_qc = False
    barcode_with_single_guide = False
    guide = ''
    cell_multiplet = False
    had_to_adjust_threshold = False
    unidentifiable_barcodes = False
    
    # NOTE TO COME BACK AND CHANGE THIS
    matches = dictionary.loc[dictionary.sorted_observable_bc_string == observed_str]

    if len(matches) > 0: # check whether the barcodes were observed at all in the table
        observed_in_dict = True
        if len(matches) == 1: # this bc combo only corresponds with one plasmid in the dictionary
            observed_barcodes_unique_plasmid = True
            barcode_passes_plasmid_qc, barcode_with_single_guide, guide = check_unambiguous_bc_combo(observed_str, matches)
        else: # the bc combo could have come from more than one plasmid
            matches_min_reads = matches.loc[matches.sufficient_reads] # filter on reads
            if len(matches_min_reads) == 1: # now try with this str
                had_to_filter_low_abundance = True
                observed_barcodes_unique_plasmid = True
                barcode_passes_plasmid_qc, barcode_with_single_guide, guide = check_unambiguous_bc_combo(observed_str, matches_min_reads)
            else: # otherwise it still could have come from multiple plasmids
                ambiguous_barcode_combination = True
    else: # if not, maybe our extraction procedure was wrong, try with a different threshold
        threshold = row['Threshold']
        if threshold > 4: # can try to come back and assign later, but these are almost certainly multiplets
            cell_multiplet = True
        else: # adjust the threshold in the most likely direction (aka 3->2, 2->3, 4->3)
            if threshold == 4:
                observed_str = row['Minus 1']
            elif threshold == 3:
                observed_str = row['Minus 1']
            elif threshold == 2:
                observed_str = row['Plus 1']

            adjusted_matches = dictionary.loc[dictionary.sorted_observable_bc_string == observed_str]

            if len(adjusted_matches) == 1: # worked, yay!
                had_to_adjust_threshold = True
                observed_in_dict = True
                observed_barcodes_unique_plasmid = True
                barcode_passes_plasmid_qc, barcode_with_single_guide, guide = check_unambiguous_bc_combo(observed_str, adjusted_matches)                    
            else: # unidentifiable combination of barcodes
                unidentifiable_barcodes = True

    return {
            'observed_str': observed_str,
            'observed_in_dict': observed_in_dict,
            'observed_barcodes_unique_plasmid': observed_barcodes_unique_plasmid,
            'had_to_filter_low_abundance': had_to_filter_low_abundance,
            'ambiguous_barcode_combination': ambiguous_barcode_combination,
            'barcode_passes_plasmid_qc': barcode_passes_plasmid_qc,
            'barcode_with_single_guide': barcode_with_single_guide,
            'OligoID': guide,
            'cell_multiplet': cell_multiplet,
            'had_to_adjust_threshold': had_to_adjust_threshold,
            'unidentifiable_barcodes': unidentifiable_barcodes
            }

# updated assignment procedure based on our new dictionary structure
# notes on the different headings from this table
# 'observed_str': sorted observed combination of barcodes in the cell used for lookup
# 'observed_in_dict': whether the bc combo was in the dictionary 
# 'observed_barcodes_unique_plasmid': whether the barcode combination could correspond to one 3xBC in the plasmid pool or many (is it consistent with multiple lines in the dictionary)
# 'had_to_filter_low_abundance': whether the barcode maps uniquely to the dict only if we filter out low abundance plasmids (<20 reads)
# 'ambiguous_barcode_combination': whether the barcode combination corresponds to multiple rows of dictionary
# 'barcode_passes_plasmid_qc': whether the barcode a) has sufficient plasmid reads and b) was actually in the array
# 'barcode_with_single_guide': whether the barcode combination corresponds to only one guide or not by plasmid sequencing
# 'OligoID': guide oligoID from the design file
# 'cell_multiplet': does this bead barcode look like it corresponds to more than one cell
# 'had_to_adjust_threshold': did I have to take +/- 1 barcode (from hypr) to get the thing to map
# 'unidentifiable_barcodes': whether the observed hypr barcodes are unidentifiable
def assign_guide_from_crop_bcs(observed_barcodes, dictionary, plots, stats, log):
    # initialize space to store the results per cell
    cell_barcodes = []
    assignment_truth_tables = []

    for barcode, row, in observed_barcodes.iterrows():
        cell_barcodes.append(barcode)
        assignment_truth_tables.append(attempt_to_assign_one_barcode(row, dictionary))

    # construct the table
    barcode_assignments = pd.DataFrame(assignment_truth_tables)
    barcode_assignments['CellBarcode'] = cell_barcodes

    # write stats
    stats.write('fraction_unique_barcode_combo\t{}\t\n'.format((barcode_assignments['observed_barcodes_unique_plasmid']).mean()))
    stats.write('fraction_single_guide\t{}\t\n'.format((barcode_assignments['OligoID'] != '').mean()))
    stats.write('fraction_real_plasmid_multiple_guides\t{}\t\n'.format(((barcode_assignments['barcode_passes_plasmid_qc']) & (barcode_assignments['OligoID'] == '')).mean()))
    stats.write('fraction_ambiguous_plasmid_given_barcode_combo\t{}\t\n'.format((barcode_assignments['ambiguous_barcode_combination']).mean()))
    stats.write('fraction_multiplets\t{}\t\n'.format((barcode_assignments['cell_multiplet']).mean()))
    stats.write('fraction_unidentifiable\t{}\t\n'.format((barcode_assignments['unidentifiable_barcodes']).mean()))
    stats.write('total_cells_assigned\t{}\t\n'.format((barcode_assignments['OligoID'] != '').sum()))

    # return the table
    return barcode_assignments

def assign_one_barcode_new(cell_barcode, observed_bcs, threshold, dictionary):
    consistent = dictionary.loc[(dictionary.sorted_bc_string==observed_bcs) | (dictionary.sorted_observable_bc_string==observed_bcs)]
    
    if len(consistent) == 1:
        assignment_type = 'single_guide'
        guide = consistent['guide'].values[0]
        barcode_combo = consistent['sorted_bc_string'].values[0]
    elif len(consistent) == 0:
        if threshold > 3:
            assignment_type = 'multiplet'
            guide = None
            barcode_combo = None
        else:
            assignment_type = 'not_in_dictionary'
            guide = None
            barcode_combo = None
    else:
        assignment_type = 'ambiguous_plasmid_origin'
        guide = ';'.join([g for g in consistent['guide'] if type(g) == str])
        barcode_combo = None

    return {'CellBarcode': cell_barcode,
            'OligoID': guide,
            'AssignmentType': assignment_type,
            'BarcodeCombination': barcode_combo}

def assign_guide_from_crop_bcs_new(observed_barcodes, dictionary, plots, stats, log):
    barcode_assignments = [] # list for collecting dicts

    for barcode, row, in observed_barcodes.iterrows():
        barcode_assignments.append(assign_one_barcode_new(barcode, row['First Choice'], row['Threshold'], dictionary))

    return pd.DataFrame(barcode_assignments) #{'OligoID': OligoIDs, 'CellBarcode': cell_barcodes})

def calculate_crop_assignment_stats(barcode_assignments, stats):
    stats.write('assigned_guides\t{}\n'.format((barcode_assignments['AssignmentType'] == 'single_guide').mean()))
    stats.write('cell_multiplets\t{}\n'.format((barcode_assignments['AssignmentType'] == 'multiplet').mean()))
    stats.write('unidentifiable_barcode_combos\t{}\n'.format((barcode_assignments['AssignmentType'] == 'not_in_dictionary').mean()))
    stats.write('ambiguous_barcode_combos\t{}\n'.format((barcode_assignments['AssignmentType'] == 'ambiguous_plasmid_origin').mean()))
    stats.write('total_cells_assigned\t{}\n'.format((barcode_assignments['AssignmentType'] == 'single_guide').sum()))


# attempts to assign a guide for each observed set of barcodes (one row in the observed barcodes file)
## DEPRECATED
def assign_guide_from_crop_bcs_old(observed_barcodes, bc_to_guide_file, plots, stats, log):
    bc_to_guide = pd.read_table(bc_to_guide_file, index_col=0)
    observed_barcodes['guide'] = observed_barcodes.apply(lambda row: get_matching_guide_from_two_dfs(row, bc_to_guide, log), axis=1)
    stats.write('Fraction assigned\t{}\n'.format((observed_barcodes['guide'] != 'No guide found').mean()))
    stats.write('Total cells assigned\t{}\n'.format((observed_barcodes['guide'] != 'No guide found').sum()))
    return observed_barcodes

def annotate_crop_guide(observed_barcodes, design_file, enhancerlist_file, promoter_file, plots, stats, log):
    design = pd.read_table(design_file)
    design = design.loc[design.subpool == 'GATA1HyPR']
    promoter = pd.read_table(promoter_file)

    # make these merges?
    observed_barcodes['element'] = observed_barcodes.apply(lambda row: design.loc[design.OligoID == row['guide'], 'target'].values[0] if row['guide'] != 'No guide found' else 'No element found', axis=1)
    observed_barcodes['TSS_symbol'] = observed_barcodes.apply(lambda row: promoter.loc[promoter.OligoID == row['guide'], 'symbol'].values[0] if row['guide'] in promoter['OligoID'].values else 'Not promoter element', axis=1)

    return observed_barcodes


# gets a standard probe name for each barcode probe (this may have to change)
def standardize_probe_names(probe_name):
    return probe_name.replace('00', '').replace('5P', '')

# valid_barcodes is list of sorted strings, counts is cells x probes (numbered)
# can optionally add guide identity too
def assign_crop_barcodes(counts, valid_barcodes, chastity, barcode_to_guide=None):
    barcodes_to_return = pd.DataFrame(index=counts.index, columns=['BarcodeStr'])

    # assign barcode strings
    for cell in counts.index:
        top_four_barcodes = counts.loc[cell].nlargest(4)
        if top_four_barcodes[2] > 0 and top_four_barcodes[2] / max(top_four_barcodes[3], 1) > chastity:
            bc_string = '-'.join(sorted(map(lambda name: standardize_probe_names(name), top_four_barcodes.index[:3].tolist())))
            if bc_string in valid_barcodes:
                barcodes_to_return.loc[cell, 'BarcodeStr'] = bc_string
            else: 
                barcodes_to_return.loc[cell, 'BarcodeStr'] = 'Passed chastity, not in bc_list'
        else:
            barcodes_to_return.loc[cell, 'BarcodeStr'] = 'Undetermined'

    # if necessary, map to their identities
    if barcode_to_guide:
        barcodes_to_return['Identity'] = barcodes_to_return.apply(lambda row: barcode_to_guide.get(row.BarcodeStr, 'Unassigned'), axis=1)

    return barcodes_to_return


# prunes a list of cbcs by removing ones that share edit distance
def remove_duplicate_cbcs(cbcs_in_cluster, led_threshold=1):
    duplicates = set()

    for i in range(len(cbcs_in_cluster)):
        for j in range(i+1, len(cbcs_in_cluster)):
            cbc1, cbc2 = cbcs_in_cluster[i], cbcs_in_cluster[j]
            if distance.levenshtein(cbc1, cbc2) <= led_threshold:
                duplicates.add(cbc2)

    return [cbc for cbc in cbcs_in_cluster if cbc not in duplicates]

# determines whether a cluster from the connected components is a singleton, multiplet, or what
# return types:
# singleton
# multiplet
# giant
def get_droplet_type(cbcs_in_cluster, small_droplet_max=3):
    if len(cbcs_in_cluster) == 1:
        return 'singleton'
    else:
        # dedup cbcs by hamming distance, returned shortened list with duplicates removed
        unique_cbcs_in_cluster = remove_duplicate_cbcs(cbcs_in_cluster)

        if len(unique_cbcs_in_cluster) == 1:
            return 'singleton'
        elif len(unique_cbcs_in_cluster) <= small_droplet_max:
            return 'multiplet'
        else:
            return 'giant'


        # # there must be a better way to do this, like what if we have two errors from the right guy...
        # if max([distance.levenshtein(cbc1, cbc2) for (cbc1, cbc2) in itertools.combinations(cbcs_in_cluster, 2)]) == 1:
        #     # droplets with more than one barcode but actually just one bead, we think
        #     return 'singleton+'
        # else:
        #     return 'multiplet'

def get_and_plot_bfp_threshold(collapsed_counts, plots, BFP_col='BFP', q=0.15):
    bfp_thresh = collapsed_counts[BFP_col].quantile(q)
    data = np.asarray(collapsed_counts[BFP_col])
    outlier_aware_hist(data, *calculate_bounds(data), cdf=False)
    plt.title('BFP Expression')
    plt.xlabel('BFP UMIs')
    plt.ylabel('Number of cells')
    plt.axvline(bfp_thresh, label='BFP threshold', c='r', ls='--')
    plt.legend()
    plt.tight_layout()
    plots.savefig()
    plt.close()
    return bfp_thresh

# filters a count matrix for just rows corresponding to singleton beads (or small no. beads)
def return_droplets_with_single_beads(counts, bc_stats, multiplets=False):
    if multiplets:
        return counts.loc[counts.index.intersection(bc_stats.loc[(bc_stats.Singleton == True) | (bc_stats.Multiplet == True)].index)]
    else:
        return counts.loc[counts.index.intersection(bc_stats.loc[bc_stats.Singleton == True].index)]

# utility function, only useful really for doing stuff in notebooks, can pass in just a count table
def return_singletons_from_counts_path(path_to_counts, multiplets=False, flat=False):
    if flat:
        counts = pd.read_table(path_to_counts).pivot(index='gene', columns='cell', values='count').fillna(0).T
    else:
        counts = pd.read_table(path_to_counts, index_col=0)
    path_to_bc_stats = glob(path.join(path.dirname(path_to_counts), 'tmp', '*.collapse_bc_stats.txt'))[0]
    bc_stats = pd.read_table(path_to_bc_stats, index_col=0)
    return return_droplets_with_single_beads(counts, bc_stats, multiplets)

# takes a list of paths to count tables (easy to use glob to get this list, probably better than just passing in directories, since we might want to do _raw or not)
# and grabs their corresponding barcode stats to get just singletons, then merges 
def join_separately_indexed_count_tables(list_of_count_tables, multiplets=False, flat=False):
    tbls = []
    for tbl in list_of_count_tables:
        tbls.append(return_singletons_from_counts_path(tbl, multiplets=multiplets, flat=flat))
    return pd.concat(tbls).fillna(0)

def rescale_matrix_to_constant_row_sums(df, use_average=True, row_sums=1000):
    average_sum = df.sum(axis=1).mean()
    if use_average:
        return df.div(df.sum(axis=1), axis=0) * average_sum
    else:
        return df.div(df.sum(axis=1), axis=0) * row_sums

# add in stuff here to compute the average correlation between probes of the same gene
# NOTE THAT FOR PLOTTING PURPOSES HERE WE ARE CLIPPING TO 0... is this right?!
def compute_and_plot_correlation(rescaled_df, plots, corr_file):
    rescaled_df_ordered = rescaled_df.reindex(sorted(rescaled_df.columns), axis=1)
    corr = rescaled_df_ordered.corr()
    fig = plt.figure(figsize=(30, 30), dpi= 80, facecolor='w', edgecolor='k')
    sns.heatmap(corr.clip(lower=0), xticklabels=corr.columns, yticklabels=corr.columns)
    plt.tight_layout()
    plt.savefig(corr_file)
    plots.savefig()
    plt.close()

    return corr

# compute the average correlation for probes within a gene and then compare that to expression data
def average_correlation_and_compare_to_expression(corr, probe_counts, plots):
    mean_probe_expression = probe_counts.mean()

    genes, probe_nums, corrs, probe_names, expressions = [], [], [], [], []

    probe_number_to_dict = {}

    for probe in corr.columns:
        gene = extract_base_probe_name(probe)
        subcorr = get_probes_for_gene(corr, gene)
        to_average = subcorr.loc[probe, [d for d in subcorr.columns if d != probe]]
        num = probe_number_to_dict.get(gene, 0)
        probe_number_to_dict[gene] = num + 1
        genes.append(gene)
        probe_nums.append(num)
        corrs.append(to_average.mean())
        probe_names.append(probe)
        expressions.append(mean_probe_expression.loc[probe])
        
    df = pd.DataFrame({'gene': genes, 'num': probe_nums, 'corr': corrs, 'probe': probe_names, 'expression': expressions})

    df['log10_expression'] = np.log10(df['expression'])

    fig, ax = plt.subplots(figsize=(25,15))
    sns.barplot(x="gene", y="corr", hue="num", data=df, ax=ax, dodge=True)
    plt.xticks(rotation=90)
    plt.ylabel('Average correlation to probes against same gene')
    plt.title('Average correlation within a gene for all probes')
    plt.tight_layout()
    plots.savefig()
    plt.close()

    df.plot.scatter('log10_expression', 'corr')
    plt.title('Per probe correlation vs mean expression')
    plt.tight_layout()
    plots.savefig()
    plt.close()

# sums across cells for either probes or genes, plots and saves the result
def sum_across_cells_and_plot(count_matrix, plots, output_file, excel_file, title='Probes'):
    marginal_counts = count_matrix.sum()

    marginal_counts = marginal_counts.loc[marginal_counts > 1]

    marginal_counts.to_csv(output_file, sep='\t')

    if excel_file:
        pd.DataFrame(marginal_counts).reset_index().to_excel(excel_file, header=False, index=False)

    marginal_counts.plot.bar()
    plt.yscale('log', nonpositive='clip')
    plt.ylabel('Counts across all cells: {}'.format(title))
    plt.xlabel(title)
    plt.title('Sum of {} across cells'.format(title))
    plt.tight_layout()
    plots.savefig()
    plt.close()

# sums across cells for either probes or genes, plots and saves the result
def average_across_cells_and_plot(count_matrix, plots, output_file, title='Probes', alpha=0.05):
    average_stats = pd.DataFrame({'mean': count_matrix.mean(), 'std': count_matrix.std(), 'sem': count_matrix.sem(), 'n': len(count_matrix), 'conf': count_matrix.sem()*norm.ppf(1-alpha/2)})

    average_stats.to_csv(output_file, sep='\t')
    # write this to excel too for jamie?

    melted = count_matrix.melt(var_name='gene')

    fig, ax = plt.subplots(figsize=(25,15))
    sns.barplot(y='value', x='gene', data=melted, ax=ax)
    plt.xticks(rotation=90)
    plt.yscale('log', nonpositive='clip')
    plt.ylabel('Average across all cells: {}'.format(title))
    plt.xlabel(title)
    plt.title('Average of {} across cells'.format(title))
    plt.tight_layout()
    plots.savefig()
    plt.close()

# returns a column of the target gene relative to the control gene(s) normalized by this ratio in negative control cells
# note that in order to be extensible this count table has to have a column called 'target' where negative control cells
# have the value 'negative_control'
def compute_normalized_expression(counts, target_gene, control_genes):
    # counts['relative'] = counts[target_gene].sum(axis=1) / counts[control_genes].sum(axis=1)
    counts['relative'] = counts[target_gene] / counts[control_genes].sum(axis=1)
    counts['normalized'] = counts['relative'] / counts.loc[counts['target'] == 'negative_control', 'relative'].mean()
    return counts['normalized']

# p16_order = ['NC-1','NC-2','NC-4','NC-5','NC-6','NC-7','NC-8','NC-9','GATA1-TSS-1','GATA1-TSS-2','HDAC6-TSS-1','HDAC6-TSS-2','e-GATA1-1','e-GATA1-2','e-HDAC6-1','e-HDAC6-2']
# p16_hue_order = ['NC', 'GATA1-TSS', 'HDAC6-TSS', 'e-GATA1', 'e-HDAC6']

def plot_p16_barplot_violinplot(counts, target_gene, plots, order=None, hue_order=None):
    # bar plot
    plt.figure(figsize=(16, 6))
    ax = sns.barplot(x='guide', y='{}_normalized'.format(target_gene), data=counts.loc[(counts.Barcode != 'Unassigned')], hue='type', dodge=False, order=order, hue_order=hue_order)
    ax.set(xticklabels=[])
    ax.set_title(target_gene)
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # violin plot with scatter
    plt.figure(figsize=(16, 6))
    ax = sns.violinplot(x='guide', y='{}_normalized'.format(target_gene), data=counts.loc[(counts.Barcode != 'Unassigned')], hue='type', dodge=False, order=order, hue_order=hue_order)
    ax = sns.stripplot(x='guide', y='{}_normalized'.format(target_gene), data=counts.loc[(counts.Barcode != 'Unassigned')], hue='type', dodge=False, jitter=True, alpha=0.3, order=order, hue_order=hue_order)
    ax.set(xticklabels=[])
    ax.set_title(target_gene)
    plt.tight_layout()
    plots.savefig()
    plt.close()


def groupby_and_report_stats(counts, target_gene, output_file, by='guide', alpha=0.05):
    counts.groupby(by)['{}_normalized'.format(target_gene)].agg([('mean', np.mean), ('std', np.std), ('n', len), ('sem', sem), ('conf_{}'.format(1-alpha), lambda value: norm.ppf(1-alpha/2) * sem(value)), ('p_val', lambda value: ttest_ind(value, counts.loc[counts.target == 'negative_control', '{}_normalized'.format(target_gene)], nan_policy='omit').pvalue)]).to_csv(output_file, sep='\t')


def compute_plot_and_save_p16_knockdown(counts, target_gene, housekeeping_genes, plots, guide_stats, element_stats, log):
    # calculate normalized column
    counts['{}_normalized'.format(target_gene)] = compute_normalized_expression(counts, target_gene, housekeeping_genes)

    # plot bar/violin plots
    plot_p16_barplot_violinplot(counts, target_gene, plots, order=sorted(counts['guide'].unique()), hue_order=sorted(counts['type'].unique()))

    # report stats for guides and elements
    groupby_and_report_stats(counts, target_gene, guide_stats, by='guide')
    groupby_and_report_stats(counts, target_gene, element_stats, by='type')


# main function that analyzes p16 data
def analyze_p16_data(collapsed_counts, housekeeping_genes, target_genes, plots, guide_stats, element_stats, log):
    # gapdh_col = collapsed_counts.filter(like='GAPDH').columns[0] # in case the GAPDH column has a new name
    # gata1_col = 'GATA1'
    # intron_cols = ['GATA1Intron', 'GAT1Int']

    # compute_plot_and_save_p16_knockdown(collapsed_counts, gata1_col, [gata1_col], gapdh_col, plots, guide_stats, element_stats, log)

    for target in target_genes:
        guide_stats_gene, element_stats_gene = guide_stats.replace('guide_stats', 'guide_stats_{}'.format(target)), element_stats.replace('element_stats', 'element_stats_{}'.format(target))
        compute_plot_and_save_p16_knockdown(collapsed_counts, target, housekeeping_genes, plots, guide_stats_gene, element_stats_gene, log)

    # # now do generically for all intron probes together
    # intron_cols_to_include = [c for c in intron_cols if c in collapsed_counts.columns]
    # log.write('Intron col: {}\n'.format('gata1-intron'))
    # guide_stats_intron, element_stats_intron = guide_stats.replace('guide_stats', 'guide_stats_gata1-intron'), element_stats.replace('element_stats', 'element_stats_gata1-intron')
    # compute_plot_and_save_p16_knockdown(collapsed_counts, 'gata1-intron', intron_cols_to_include, gapdh_col, plots, guide_stats_intron, element_stats_intron, log)
    
    # # now do HDAC6 + intron
    # guide_stats_hdac6, element_stats_hdac6 = guide_stats.replace('guide_stats', 'guide_stats_hdac6'), element_stats.replace('element_stats', 'element_stats_hdac6')
    # compute_plot_and_save_p16_knockdown(collapsed_counts, 'HDAC6', ['HDAC6NM_001321225'], gapdh_col, plots, guide_stats_hdac6, element_stats_hdac6, log)

    # guide_stats_hdac6_int, element_stats_hdac6_int = guide_stats.replace('guide_stats', 'guide_stats_hdac6_int'), element_stats.replace('element_stats', 'element_stats_hdac6_int')
    # compute_plot_and_save_p16_knockdown(collapsed_counts, 'HDAC6-intron', ['HDAC6Intron4'], gapdh_col, plots, guide_stats_hdac6_int, element_stats_hdac6_int, log)

def plot_crop_dictionary(dictionary, plots, reads_thresh=None, fraction_thresh=None):
    # plt.scatter(np.log10(dictionary['reads']), dictionary['fraction'], alpha=0.1)
    plt.scatter(np.log10(dictionary.loc[dictionary.guide_mapped, 'reads']), dictionary.loc[dictionary.guide_mapped, 'fraction'], alpha=0.1, c='blue', label='assigned')
    plt.scatter(np.log10(dictionary.loc[~dictionary.guide_mapped, 'reads']), dictionary.loc[~dictionary.guide_mapped, 'fraction'], alpha=0.1, c='red', label='unassigned')
    plt.xlabel('Reads per barcode triplet (log10)')
    plt.ylabel('Fraction reads supporting top guide')
    if reads_thresh:
        plt.axvline(np.log10(reads_thresh), c='k', ls='--', label='Min reads: {}'.format(reads_thresh))
    if fraction_thresh:
        plt.axhline(fraction_thresh, c='k', ls='--', label='Fraction: {:.2f}'.format(fraction_thresh))
    plt.legend()
    plt.tight_layout()
    plots.savefig()
    plt.close()

# return a list of genes that are "expressed" from a probe or gene count matrix
# only returns genes that are present over a certain expression level
# can optionally remove barcode probes or blacklist specific probes/genes
def get_expressed_genes(count_table, remove_barcodes=False, blacklist=None, whitelist=None, min_expression=1):
    all_genes = count_table.columns.tolist()

    count_table_rescaled = rescale_matrix_to_constant_row_sums(count_table, use_average=True)

    filtered_genes = [g for g in all_genes if count_table_rescaled[g].mean() > min_expression]

    if remove_barcodes:
        filtered_genes = [g for g in filtered_genes if 'bar' not in g]

    if whitelist:
        filtered_genes = set(filtered_genes + whitelist)

    if blacklist:
        filtered_genes = [g for g in filtered_genes if g not in blacklist]

    return filtered_genes

# uses zero truncated poisson dist (above) to get mle on cell lambda
# this is specifically for p16 data, but works the same way as the other ones, above
def calculate_cell_lambda_from_p16_counts(counts, umi_thresh=10):
    barcode_counts = counts.filter(like='bar')

    cells_per_droplet = (barcode_counts > umi_thresh).sum(axis=1)

    n_cells = (cells_per_droplet > 0).sum() # just look at cells with at least one guide

    counts_per_droplet_type = cells_per_droplet.value_counts().reset_index()

    x_bar = (counts_per_droplet_type['index'] * counts_per_droplet_type[0]).sum() / n_cells

    optim_res = minimize_scalar(lambda x: zero_trunc_poisson_objective(x, x_bar))

    if optim_res.success:
        return optim_res.x
    else:
        return -1

# normalize a count table to the housekeeping genes per row
def normalize_counts_to_housekeeping(cts, hk_list, scaling_factor=1):
    return cts.div(cts[hk_list].sum(axis=1), axis=0)*scaling_factor

def remove_nm_number(base_probe_name):
    if 'NM' in base_probe_name:
        return base_probe_name.split('_')[0][:-2]
    else:
        return base_probe_name

# turn count matrix (probes) into count matrix (genes), potentially excluding bad probes
def sum_probes_with_exclusion_list(probe_counts, exclusion_list=[]):
    gene_counts = pd.DataFrame(index=probe_counts.index)
    
    genes = sorted(list({extract_base_probe_name(probe) for probe in probe_counts.columns}))
    
    for gene in genes:
        # need to add a check in here if I want to use this in analyze_count for the barcode probes
        gene_counts[gene] = get_probes_for_gene(probe_counts, gene, exclusion_list).sum(axis=1)
        
    return gene_counts

# identify probes that are far from the median probe, plot them, and return a list
def find_problematic_probes(probes, plots, stats, stats_suffix, problematic_probes_file, threshold=1):
    to_plot = pd.DataFrame(columns=['gene', 'log2(fold change) relative to median'])

    # make sure this works for singleton probes, barcode probes
    for g in { extract_base_probe_name(probe) for probe in probes.columns }:
        avg = get_probes_for_gene(probes, g).mean()
        avg /= np.median(avg.values)
        df = pd.DataFrame({'Gene': g, 'log2(fold change) relative to median': np.log2(avg)})
        to_plot = to_plot.append(df)

    # plot
    fig, ax = plt.subplots()
    gene_order = sorted(to_plot.Gene.unique())
    
    sns.stripplot(data=to_plot, x='Gene', y='log2(fold change) relative to median', ax=ax, size=15, color='k', jitter=False, order=gene_order)
    ax.axhline(threshold, c='k', ls='--', lw=3)
    ax.axhline(-threshold, c='k', ls='--', lw=3)

    ax.set_xlabel('')
    plt.xticks(rotation=90)

    # max_val, min_val = to_plot['log2(fold change) relative to median'].dropna().values.max(), to_plot['log2(fold change) relative to median'].dropna().values.min()
    # plt.yticks(range(int(min_val), 2+int(max_val)))

    plt.tight_layout()
    plots.savefig()
    plt.close()

    # compute stats
    good_probes = to_plot['log2(fold change) relative to median'].between(-threshold, threshold, inclusive=True)
    fraction_good_probes = good_probes.mean()
    stats.write('Fraction_good_probes{}\t{}\n'.format(stats_suffix, fraction_good_probes))
    stats.write('Probe_threshold{}\t{}\n'.format(stats_suffix, threshold))

    # write problematic probes to file
    to_plot.loc[~good_probes].to_csv(problematic_probes_file, sep='\t', header=True, index=True)
