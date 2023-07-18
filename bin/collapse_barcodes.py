#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.sparse import csc_matrix, save_npz
from scipy.special import comb
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import argparse

import mip_tools


# read in groups file, analyze the overlap between barcodes, and write a new whitelist to re-extract
def collapse_beads_and_write_whitelist(groups_file, whitelist_file, filtered_whitelist_file, barcode_stats, overlap_mat, plots, log, stats, thresh, histogram_file, clustermap_png):
    # read in file
    groups = pd.read_table(groups_file, usecols=['read_id', 'umi', 'final_umi', 'contig', 'unique_id'])
    
    
    # annotate groups df with additional columns
    groups['CBC'] = groups.apply(lambda row: mip_tools.extract_cbc(row['read_id']), axis=1)
    groups['pathological'] = groups.apply(lambda row: mip_tools.is_pathological(row['CBC']), axis=1)
    groups['UMI_Transcript'] = groups['final_umi'] + '_' + groups['contig']

    # filter groups before proceeding to calculate overlap
    ### QUESTION -- definitely need to filter out the G barcode here, but what about the other pathological ones?
    filtered = groups.loc[~groups.pathological]
    print("Filtering done")

    mip_tools.make_knee_plot_from_groups(filtered, plots)

    # write how many BCs dropped for being pathological
    path_cbcs = len(groups.loc[groups.pathological, 'CBC'].unique())
    stats.write('Pathological CBCs\t{}\n'.format(path_cbcs))

    # write some stats
    stats.write('Reads (groups file)\t{}\n'.format(len(groups)))
    stats.write('Reads (filtered for pathological BCs)\t{}\n'.format(len(filtered)))
    # additional states: total reads, total CBCs, path CBCs, etc. etc. (as much as possible, make sure to specify where)

    # build adjacency matrix ( and get stats)
    adjacency_matrix, cbcs, n_beads, sets = mip_tools.build_adjacency_matrix(filtered)
    avg_overlap_stat = adjacency_matrix.mean() / (1 - 1/adjacency_matrix.shape[0]) # correction for the diagonals
    print("Built adj matrix")
    # save adjacency matrix
    save_npz(overlap_mat, csc_matrix(adjacency_matrix))

    # binarize adjacency matrix
    adjacency_matrix_filtered = mip_tools.binarize_adjacency_matrix(adjacency_matrix, thresh)

    # do connected components
    n_clusters, labels = mip_tools.do_connected_components(adjacency_matrix_filtered)
    print("Find connected components")
    # plot adjacency matrix (now taking into account the cluster assignments)
    mip_tools.plot_adjacency_matrix(adjacency_matrix, labels, plots, thresh, log, histogram_file, clustermap_png)

    # write and plot results
    stats.write('Number of clusters\t{}\n'.format(n_clusters))
    stats.write('Avg overlap stat (log)\t{}\n'.format(np.log(avg_overlap_stat)))
    # mip_tools.plot_barcodes_per_droplet(labels, plots)

    # initialize the barcode stats df
    # bc_stats = pd.DataFrame(index=cbcs, columns=['Pathological', 'Singleton', 'Singleton+', 'Multiplet', 'Representative', 'n_cbcs', 'n_beads'])
    bc_stats = pd.DataFrame(index=cbcs, columns=['Pathological', 'Singleton', 'Multiplet', 'Giant', 'ID', 'n_cbcs', 'n_beads'])
    # instead of this representative barcode thing (which we should coordinate with whitelist, probably, instead just include a unique ID per cluster, will make stats below much easier)

    # first add pathological ones (with stuff here to show that they are not included in the clusters, aka the -1)
    for path_bc in groups.loc[groups.pathological, 'CBC'].unique():
        path, sing, mult, giant = True, False, False, False
        bc_stats.loc[path_bc] = [path, sing, mult, giant, -1, 1, 1]

    # iterate through the clusters to build the df (and whitelist)
    whitelist = []
    filtered_whitelist = []

    for idx in range(n_clusters):
        cbcs_in_cluster = cbcs[labels == idx].tolist()

        bc_type = mip_tools.get_droplet_type(cbcs_in_cluster)

        whitelist.append(cbcs_in_cluster)

        if bc_type in ['singleton', 'multiplet']:
            filtered_whitelist.append(cbcs_in_cluster.copy())

        if bc_type == 'singleton':
            path, sing, mult, giant = False, True, False, False

        elif bc_type == 'multiplet':
            path, sing, mult, giant = False, False, True, False

        elif bc_type == 'giant':
            path, sing, mult, giant = False, False, False, True

        else: # so far nothing
            path, sing, mult, giant = False, False, False, False

        n_cbcs = len(cbcs_in_cluster)
        n_beads = len(mip_tools.remove_duplicate_cbcs(cbcs_in_cluster))

        for bcx in cbcs_in_cluster:
            bc_stats.loc[bcx] = [path, sing, mult, giant, idx, n_cbcs, n_beads]


                # bc = cbcs_in_cluster[0]
            # bc_stats.loc[bc] = [False, True, False, False, bc, 1, 1]

        # elif bc_type == 'singleton+':
        #     bc1 = cbcs_in_cluster[0]
        #     for bcx in cbcs_in_cluster:
        #         bc_stats.loc[bcx] = [False, True, True, False, bc1, len(cbcs_in_cluster), 1]

        # else: # multiplet
        #     bc1 = cbcs_in_cluster[0]
        #     for bcx in cbcs_in_cluster:
        #         bc_stats.loc[bcx] = [False, False, False, True, bc1, len(cbcs_in_cluster), len(cbcs_in_cluster)]

    # plot stats about the barcodes
    # mip_tools.plot_barcode_stats(bc_stats, plots)

    # write stats
    stats.write('Number of singleton clusters\t{}\n'.format(len(bc_stats.loc[bc_stats['Singleton']==True, 'ID'].unique())))
    stats.write('Number of singleton barcodes\t{}\n'.format(sum(bc_stats['Singleton'])))
    stats.write('Number of multiplet clusters\t{}\n'.format(len(bc_stats.loc[bc_stats['Multiplet']==True, 'ID'].unique())))
    stats.write('Number of multiplet barcodes\t{}\n'.format(sum(bc_stats['Multiplet'])))
    stats.write('Number of giant clusters\t{}\n'.format(len(bc_stats.loc[bc_stats['Giant']==True, 'ID'].unique())))
    stats.write('Number of giant barcodes\t{}\n'.format(sum(bc_stats['Giant'])))
    stats.write('Fraction singleton barcodes\t{}\n'.format(np.mean(bc_stats['Singleton'])))
    # stats.write('Number of singleton clusters\t{}\n'.format(len(bc_stats.loc[(bc_stats['Singleton']==True) & ~(bc_stats['Singleton+']==True)]) + len(bc_stats.loc[bc_stats['Singleton+']==True, 'Representative'].unique())))
    # stats.write('Number of singleton barcodes\t{}\n'.format(sum(bc_stats['Singleton'])))
    # stats.write('Number of multiplet clusters\t{}\n'.format(len(bc_stats.loc[bc_stats['Multiplet']==True, 'Representative'].unique())))
    # stats.write('Number of multiplet barcodes\t{}\n'.format(sum(bc_stats['Multiplet'])))
    # stats.write('Percent singleton barcodes\t{}\n'.format())
    stats.write('Percent barcode pairs above threshold\t{:.3f}\n'.format(100 * np.count_nonzero(adjacency_matrix_filtered)/2 / comb(len(cbcs), 2)))

    est_lambda = mip_tools.estimate_lambda(bc_stats)

    stats.write('Estimated bead lambda\t{:.3f}\n'.format(est_lambda))

    # write whitelists
    mip_tools.write_whitelist(whitelist, whitelist_file, sets)
    mip_tools.write_whitelist(filtered_whitelist, filtered_whitelist_file, sets)

    # write stats df
    bc_stats.to_csv(barcode_stats, sep='\t')



    # # initialize data structures to hold this info
    # whitelist = [] # list of bc clusters: [[BC1], [BC2,BC3], ...]
    # singletons = [] # list of bcs that are singletons (if two differ by 1 led, we just choose the first one in WL)
    # multiplets = [] # list of bcs are are multiplets (always choose representative of cluster, same as whitelist)

    # # iterate through all clusters, writing both the whitelist and the singlet/multiplet files
    # for idx in range(n_clusters):
    #     cbcs_in_cluster = cbcs[labels == idx]
    #     if len(cbcs_in_cluster) == 1:
    #         cbc = cbcs_in_cluster[0]
    #         whitelist.append([cbc])
    #         singletons.append(cbc)
    #     else:
    #         # get the list of all the cbcs in the cluster
    #         cbcs_in_cluster = cbcs_in_cluster.tolist()

    #         # calculate the maximum edit distance between all the pairs of cbcs
    #         # if it's 1 we can treat it as actually a singlet
    #         ### NOTE: PROBABLY A BETTER WAY TO DO THIS IS WHETHER WE CAN FIND A 1 LED PATH ALL THE WAY THROUGH THE 
    #         ### GRAPH, BUT THIS SHOULD WORK FOR NOW
    #         min_edit_distance = min([distance.levenshtein(cbc1, cbc2) for (cbc1, cbc2) in itertools.combinations(cbcs_in_cluster, 2)])

    #         # # grab representative cbc and remainder
    #         # cbc = cbcs_in_cluster.pop(0)

    #         # # add to lists
    #         # if min_edit_distance == 1:
    #         #     singletons.append(cbc)
    #         # else:
    #         #     multiplets.append(cbc)

    #         # whitelist.append([cbc]+cbcs_in_cluster)

    #         # making a change here to actually write all the barcodes to the multiplets and singletons file (so that we can color in P16)
    #         # presumably we should write two different lists, one for the collapsed index and one for the non-, but we can add this later

    #         # add to lists
    #         if min_edit_distance == 1:
    #             for c in cbcs_in_cluster:
    #                 singletons.append(c)
    #         else:
    #             for c in cbcs_in_cluster:
    #                 multiplets.append(c)

    #         whitelist.append(cbcs_in_cluster)
            

    # # write stats on barcode distributions
    # multiplet_bcs = [multi_lst for multi_lst in whitelist if len(multi_lst) > 1] # this is everything with more than 1 bc, not necessarily more than 1 bead though... 
    # stats.write('Number of singleton barcodes\t{}\n'.format(len(singletons)))
    # stats.write('Number of droplets with multiplet barcodes\t{}\n'.format(len(multiplet_bcs)))
    # stats.write('Number of droplets with multiplet barcodes but one bead\t{}\n'.format(len(multiplet_bcs)-len(multiplets)))
    # stats.write('Number of barcodes in multiplet droplets\t{}\n'.format(sum([len(multi) for multi in multiplet_bcs])))
    # # still some stuff to clean up around here... although the numbers to work out

    # # write out files
    # # with open(output_prefix+'_singleton_bcs.txt', 'w') as singletons_file, open(output_prefix+'_multiplet_bcs.txt', 'w') as multiplets_file:
    # ### CHANGE THIS TO INSTEAD WRITE A DF WITH ALL THE INFO FOR THE BCS, THEN REPORT STATS
    # with open(singletons_path, 'w') as singletons_file, open(multiplets_path, 'w') as multiplets_file:
    #     for bc in singletons:
    #         singletons_file.write('{}\n'.format(bc))
    #     for bc in multiplets:
    #         multiplets_file.write('{}\n'.format(bc))
    
    # mip_tools.write_whitelist(whitelist, whitelist_file, sets)

    log.write('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze a groups file and write a whitelist, taking into account collapsing overlapping barcodes')
    
    parser.add_argument("-g", "--groups", dest="groups_file", type=str,  help="Groups file")
    parser.add_argument("-w", "--whitelist", dest="whitelist_file", type=str, help="Path to output whitelist")
    parser.add_argument("-f", "--filtered_whitelist", dest="filtered_whitelist_file", type=str, help="Path to output filtered whitelist")
    # parser.add_argument("-o", "--output", dest="output_prefix", type=str, help="Path prefix to output files, for singletons and multiplet bcs")
    parser.add_argument("-b", "--barcode_stats", dest="barcode_stats", type=str, help="Path to output whitelist")
    parser.add_argument("-o", "--overlap_mat", dest="overlap_mat", type=str, help="Path to output whitelist")
    # parser.add_argument("-m", "--multiplets", dest="multiplets", type=str, help="Path to write multiplet barcodes")
    # parser.add_argument("-z", "--singletons", dest="singletons", type=str, help="Path to write singletons barcodes")
    parser.add_argument("-p", "--plots", dest="plots_file", type=str,  help="Where to output plots")
    parser.add_argument("-l", "--log", dest="log_file", type=str,  help="Where to output log")
    parser.add_argument("-s", "--stats", dest="stats_file", type=str,  help="Where to output stats")
    parser.add_argument("-t", "--thresh", dest="thresh", type=float, help="Threshold for determining overlap")
    parser.add_argument("-c", "--histogram", dest="histogram_file", type=str,  help="Where to output histogram png")
    parser.add_argument("-m", "--clustermap", dest="clustermap_png", type=str,  help="Where to output clustermap png")

    args = parser.parse_args()
    print("HELLO WORLD")
    with PdfPages(args.plots_file) as plots, open(args.stats_file, 'w') as stats, open(args.log_file, 'w') as log:
        print("COLLAPSE BARCODE IS RUNNING")
        collapse_beads_and_write_whitelist(args.groups_file, args.whitelist_file, args.filtered_whitelist_file, args.barcode_stats, args.overlap_mat, plots, log, stats, args.thresh, args.histogram_file, args.clustermap_png)

