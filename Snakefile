import pandas as pd
import numpy as np
from os import path

# monkey patch the snakemake SequenceFormatter for naming submitted jobs
from snakemake.utils import SequenceFormatter
from snakemake.io import Wildcards

original_formatter = SequenceFormatter.format_field

def new_format_field(self, value, format_spec):
    if isinstance(value, Wildcards):
        return ".".join("{}={}".format(name, value) # simply change the ',' to '.'
                        for name, value in
                        sorted(value.items(), key=lambda item: item[0]))
    if isinstance(value, (list, tuple, set, frozenset)):
        return self.separator.join(self.format_element(v, format_spec)
                                   for v in value)
    else:
        return self.format_element(value, format_spec)

SequenceFormatter.format_field = new_format_field

# defining global variables
mip_directory = '/seq/lincRNA/Projects/MIP/mip-pipeline/bin'

def parse_sample_sheet(path_to_sample_sheet):
    # add in column for downsampling, if necessary
    sample_sheet = pd.read_table(path_to_sample_sheet)
    if 'ReadsToUse' not in sample_sheet.columns:
        sample_sheet['ReadsToUse'] = False
    return sample_sheet

# getting arguments passed into config
mismatches = config['mismatches']
umi = ''.join(['N' for _ in range(config['umi'])])
method = config['method']
hamming = config['hamming']
readlen = config['readlen']
cbc_len = config['cbc']
cbc = ''.join(['C' for _ in range(cbc_len)])
thresh = config.get('thresh', 0.005) # 0.005 what about 0.01?

# get the sample sheet to determine what we actually need to run on
sample_sheet = parse_sample_sheet(config['sample_sheet'])

# converts a sample/folder name (which is NAME-ID) and just returns ID
def sample_name_to_id(sample):
    return int(sample.split('-')[-1])

# read from sample sheet and define target for all samples
rule all:
    input: "summary.txt"

# figure out whether we need to produce bulk or emul outputs, or both
def collate_helper(wildcards):
    n_bulk_samples = (sample_sheet.Read2 == 'Bulk').sum()
    targets = []
    if n_bulk_samples:
        targets.append('bulk_count_summary.txt')
    if n_bulk_samples < len(sample_sheet):
        targets.append('emul_count_summary.txt')
    return targets

# collate the information from the different summaries into one document
rule collate:
    input:
        collate_helper
    output:
        "summary.txt"
    params:
        collate=path.join(mip_directory, 'collate_hypr.py')
    shell:
        # WRITE THIS
        'python {params.collate} -o {output} -i {input}'        

# glue stats horizontally (i.e. by experiment)
rule bulk_summary:
    input:
        counts=["{name}-{id}/{name}-{id}_bulk.counts.txt".format(name=row.Name, id=row.ID) for row in sample_sheet.loc[sample_sheet.Read2 == 'Bulk'].itertuples()],
        stats=["{name}-{id}/{name}-{id}_bulk.stats.txt".format(name=row.Name, id=row.ID) for row in sample_sheet.loc[sample_sheet.Read2 == 'Bulk'].itertuples()]
        # make this require other stuff, like the other outputs of analyze group or whatever...
        # expand(["{nameid}/{nameid}_bulk.counts.txt", "{nameid}/{nameid}_bulk.stats.txt"], nameid=['{}-{}'.format(row.Name, row.ID) for row in sample_sheet.loc[sample_sheet.BarcodeRead == 'Bulk'].itertuples()])
    output:
        counts_summary='bulk_count_summary.txt',
        stats_summary='bulk_stats_summary.txt'
    params:
        bulk_summarize=path.join(mip_directory, 'bulk_summarize.py')
    shell:
        "python {params.bulk_summarize} --output_stats {output.stats_summary} --input_stats {input.stats} \
            --output_counts {output.counts_summary} --input_counts {input.counts}"

rule emul_summary:
    input:
        counts=["{name}-{id}/{name}-{id}_collapsed_emul.sum_genes.txt".format(name=row.Name, id=row.ID) for row in sample_sheet.loc[sample_sheet.Read2 != 'Bulk'].itertuples()],
        stats=["{name}-{id}/{name}-{id}.emul_stats.txt".format(name=row.Name, id=row.ID) for row in sample_sheet.loc[sample_sheet.Read2 != 'Bulk'].itertuples()]
    output:
        counts_summary='emul_count_summary.txt',
        stats_summary='emul_stats_summary.txt'
    params:
        emul_summarize=path.join(mip_directory, 'emul_summarize.py')
    shell:
        "python {params.emul_summarize} --output_stats {output.stats_summary} --input_stats {input.stats} \
            --output_counts {output.counts_summary} --input_counts {input.counts}"

# bulk specific analysis
rule bulk_output_plots:
    input:
        # groups="{sample}/tmp/{sample}_bulk.groups.tsv",
        groups="{sample}/tmp/{sample}_bulk.groups.tsv.gz",
        bwtlog="{sample}/log/{sample}_bulk.bowtie.log"
    output:
        plot="{sample}/{sample}.group_plot.pdf",
        log="{sample}/log/{sample}.group.plots.log",
        stats="{sample}/{sample}_bulk.stats.txt"
    shell:
        "python {mip_directory}/bulk_output_plots.py -g {input.groups} \
            -p {output.plot} -l {output.log} -s {output.stats} -b {input.bwtlog}"

rule parse_bulk_count:
    input:
        raw="{sample}/{sample}_bulk_raw.counts.txt"
    output:
        counts="{sample}/{sample}_bulk.counts.txt",
        plot="{sample}/{sample}.probe_plot.pdf",
        probes="{sample}/{sample}_bulk.probe_counts.txt"
    shell:
        "{mip_directory}/parse_count.py -r {input.raw} -c {output.counts} -p {output.probes} -o {output.plot}"

# emul specific analyses (glue vertically, as in combine stats from a bunch of separate places into one)
rule combine_emul_stats:
    input:
        dup_stats_un="{sample}/stats/{sample}_uncollapsed.duplication_stats.txt",
        dup_stats="{sample}/stats/{sample}_collapsed.duplication_stats.txt",
        dup_stats_filt="{sample}/stats/{sample}_collapsed_filtered.duplication_stats.txt",
        count_stats_un="{sample}/stats/{sample}_uncollapsed.out_stats.txt",
        count_stats="{sample}/stats/{sample}_collapsed.out_stats.txt",
        count_stats_filt="{sample}/stats/{sample}_collapsed_filtered.out_stats.txt",
        collapse_stats="{sample}/stats/{sample}.collapse_barcodes.stats.txt",
        collapsed_reads="{sample}/stats/{sample}_emul_extracted_collapsed.stats.txt",
        uncollapsed_reads="{sample}/stats/{sample}_emul_extracted_uncollapsed.stats.txt",
        collapsed_filt_reads="{sample}/stats/{sample}_emul_extracted_collapsed_filtered.stats.txt",
        filtered_reads="{sample}/stats/{sample}_filtered.stats.txt",
        input_reads="{sample}/stats/{sample}_downsampled.stats.txt",
        starting_cells="{sample}/stats/{sample}_starting_cells.stats.txt"
    output:
        "{sample}/{sample}.emul_stats.txt"
    shell:
        "python {mip_directory}/combine_stats.py -i {input.dup_stats} {input.dup_stats_un} {input.dup_stats_filt} {input.count_stats} {input.count_stats_un} {input.count_stats_filt} \
            {input.collapse_stats} {input.collapsed_reads} {input.uncollapsed_reads} {input.collapsed_filt_reads} {input.filtered_reads} {input.input_reads} \
            {input.starting_cells} -o {output} -e .emul"

rule analyze_count:
    input:
        flat="{sample}/{sample}_{collapsed}_emul_flat.counts.txt"
    output:
        stats="{sample}/stats/{sample}_{collapsed}.out_stats.txt",
        plots="{sample}/{sample}_{collapsed}.out_plots.pdf",
        log="{sample}/log/{sample}_{collapsed}.analyze_count.log",

        cells_by_probes="{sample}/{sample}_{collapsed}_emul.count_matrix_probes.txt",
        cells_by_probes_rescaled="{sample}/{sample}_{collapsed}_emul.count_matrix_probes_rescaled.txt",
        cells_by_genes="{sample}/{sample}_{collapsed}_emul.count_matrix.txt",
        cells_by_genes_rescaled="{sample}/{sample}_{collapsed}_emul.count_matrix_rescaled.txt",
        sum_probes_across_cells="{sample}/{sample}_{collapsed}_emul.sum_probes.txt",
        sum_genes_across_cells="{sample}/{sample}_{collapsed}_emul.sum_genes.txt",
        average_probes_across_cells="{sample}/{sample}_{collapsed}_emul.average_probes.txt",
        average_genes_across_cells="{sample}/{sample}_{collapsed}_emul.average_genes.txt",

        correlation="{sample}/{sample}_{collapsed}.probe_correlation.png",
        probes_zscores="{sample}/{sample}_{collapsed}_emul.probe_zscores.txt",
        genes_zscores="{sample}/{sample}_{collapsed}_emul.gene_zscores.txt",
        problematic_probes="{sample}/{sample}_{collapsed}_emul.problematic_probes.txt",
        
        excel_counts="outputs/probecountspercell{sample}_{collapsed}.xlsx",
        excel_probes="outputs/countsperprobe{sample}_{collapsed}.xlsx",
        excel_genes="outputs/countspergene{sample}_{collapsed}.xlsx",
        excel_zscores="outputs/zscorespercell{sample}_{collapsed}.xlsx",
    params:
        # collapsed=lambda wildcards: wildcards.collapsed == 'collapsed',
        collapsed=lambda wildcards: wildcards.collapsed,
        bc_stats=lambda wildcards: "{}/tmp/{}.collapse_bc_stats.txt".format(wildcards.sample, wildcards.sample) if '_collapsed' in wildcards.collapsed else 'None'
    shell:
         "python {mip_directory}/analyze_count.py --flat {input.flat} \
            --cell_x_gene {output.cells_by_genes} --cell_x_gene_rescaled {output.cells_by_genes_rescaled} --cell_x_probe {output.cells_by_probes} --cell_x_probe_rescaled {output.cells_by_probes_rescaled} \
            --sum_probes {output.sum_probes_across_cells} --sum_genes {output.sum_genes_across_cells} --average_probes {output.average_probes_across_cells} --average_genes {output.average_genes_across_cells} \
            --plots {output.plots} --stats {output.stats} --log {output.log} --collapsed {params.collapsed} --bc_stats {params.bc_stats} \
            --gene_z {output.genes_zscores} --probe_z {output.probes_zscores} --corr {output.correlation} --problematic_probes {output.problematic_probes} \
            --excelcounts {output.excel_counts} --excelprobes {output.excel_probes} --excelgenes {output.excel_genes} --excelzscores {output.excel_zscores}"

# analyzes grouped emul HCR data for info on bead doublets, duplication rate, etc.
rule analyze_emul_duplication:
    input:
        # groups='{sample}/tmp/{sample}_{collapsed}.groups.tsv',
        groups_zip="{sample}/tmp/{sample}_{collapsed}.groups.tsv.gz", # this ensures that we make the gzipped version
        bwtlog="{sample}/log/{sample}_{collapsed}.bowtie.log"
        # wl='{sample}/tmp/{sample}_whitelist.txt'
    output:
        log="{sample}/log/{sample}_{collapsed}.analyze_duplication.log",
        stats="{sample}/stats/{sample}_{collapsed}.duplication_stats.txt",
        plot="{sample}/{sample}_{collapsed}.duplication_stats.pdf"
    params:
        # collapsed=lambda wildcards: wildcards.collapsed == 'collapsed'
        collapsed=lambda wildcards: wildcards.collapsed
    shell:
        "python {mip_directory}/analyze_duplication.py -g {input.groups_zip} -l {output.log} \
            -s {output.stats} -b {input.bwtlog} -p {output.plot} -u {params.collapsed}"

# count rules (both emulsion and bulk)
rule bulk_count:
    input:
        "{sample}/tmp/{sample}_bulk.extracted.trim.sorted.grouped.bam"
    output:
        "{sample}/{sample}_bulk_raw.counts.txt"
    shell:
        "samtools view {input} | cut -f 3,16 | uniq | cut -f 1 | uniq -c | sort -nr > {output}"
    
rule emul_count:
    input:
        bam="{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam", #"{sample}/tmp/{sample}_collapsed.extracted.trim.sorted.bam",
        index="{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam.bai"
    output:
        counts="{sample}/{sample}_{collapsed}_emul_flat.counts.txt",
        log="{sample}/log/{sample}_{collapsed}.count.log",
    shell:
        "umi_tools count --per-gene --per-contig --per-cell --method={method} \
                   -I {input.bam} -S {output.counts} -L {output.log}"

# rule for figuring out what barcodes to include (before writing new whitelist and rerunning the pipeline)
rule collapse_barcodes:
    input:
        # groups="{sample}/tmp/{sample}_uncollapsed.groups.tsv"
        groups="{sample}/tmp/{sample}_uncollapsed.groups.tsv.gz"
    output:
        whitelist='{sample}/tmp/{sample}_collapsed_whitelist.txt',
        filtered_whitelist='{sample}/tmp/{sample}_collapsed_filtered_whitelist.txt',
        barcode_stats='{sample}/tmp/{sample}.collapse_bc_stats.txt',
        # singletons='{sample}/tmp/{sample}.collapse_singleton_bcs.txt',
        plots="{sample}/{sample}.collapse_barcodes.plots.pdf",
        log="{sample}/log/{sample}.collapse_barcodes.log.txt",
        stats="{sample}/stats/{sample}.collapse_barcodes.stats.txt",
        overlap_matrix=temp('{sample}/tmp/{sample}.collapse_adjacency_mat.npz'),
        histogram='{sample}/{sample}.collapse_barcodes_histogram.png'
    params:
        thresh=thresh,
        clustermap=lambda wildcards: '{}/{}.collapse_barcodes_clustermap_raw.png'.format(wildcards.sample, wildcards.sample)
    shell:
        "python {mip_directory}/collapse_barcodes.py -g {input.groups} -w {output.whitelist} -f {output.filtered_whitelist} -b {output.barcode_stats} \
            -o {output.overlap_matrix} -p {output.plots} -l {output.log} -s {output.stats} -t {params.thresh} -c {output.histogram} -m {params.clustermap}"
        
# common rules (most of the rules are generic for bulk and emulsion)
# only note is that whether or not to include --per-cell for grouping depends on type

# this gets rid of the bulky groups file when it is no longer needed (it is temp) and stores as gzip
# this rule doesn't actually remove the tsv (since we are using rediction), but the temp in the group
# rule should delete it when it is no longer useful
# but I think we can actually just use the gziped version everywhere, because pandas is smart

# 200106: also making gziped groups temp. this is still the vast majority of our storage, and I don't think 
# it's actually necessary to hold onto it (especially since I'm not really changing the things that use it anymore)
rule gzip_group:
    input:
        "{sample}/tmp/{sample}_{collapsed}.groups.tsv"
    output:
        # temp("{sample}/tmp/{sample}_{collapsed}.groups.tsv.gz")
        "{sample}/tmp/{sample}_{collapsed}.groups.tsv.gz"
    shell:
        "gzip < {input} > {output}"


# group reads (depends on emulsion vs bulk)
rule group:
    input:
        bam="{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam",
        index="{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam.bai"
    output:
        groupedbam=temp("{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.grouped.bam"),
        log="{sample}/log/{sample}_{collapsed}.group.log",
        groups=temp("{sample}/tmp/{sample}_{collapsed}.groups.tsv")
    params:
        percell=lambda wildcards: '' if sample_sheet.loc[sample_sheet.ID == sample_name_to_id(wildcards.sample), 'Read2'].values[0] == 'Bulk' else '--per-cell'
    shell:
        "umi_tools group {params.percell} --output-bam --group-out={output.groups} \
            -I {input.bam} -S {output.groupedbam} -L {output.log} \
            --method={method} --edit-distance-threshold={hamming}"

# build bam index
# Note, we don't actually use {output} in the shell command, but Snakemake needs to know that it creates this file
rule index_bam:
    input:
        "{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam"
    output:
        temp("{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam.bai")
    shell:
        "samtools index {input}"

# sort bam
rule sort_bam:
    input:
        "{sample}/tmp/{sample}_{collapsed}.extracted.trim.bam"
    output:
        temp("{sample}/tmp/{sample}_{collapsed}.extracted.trim.sorted.bam")
    shell:
        "samtools sort {input} -o {output}"

# compress sam file
rule sam_to_bam:
    input:
        "{sample}/tmp/{sample}_{collapsed}.extracted.trim.sam"
    output:
        temp("{sample}/tmp/{sample}_{collapsed}.extracted.trim.bam")
    shell:
        "samtools view -b {input} -o {output}"

# finds the proper bowtie index for each sample
def map_reads_helper(wildcards):
    sample = wildcards.sample
    index_name = sample_sheet.loc[sample_sheet.ID == sample_name_to_id(sample), 'ProbeName'].values[0]
    return {'trim': '{sample}/tmp/{sample}_{collapsed}_R1_001.extracted.trim.fastq.gz'.format(sample=sample, collapsed=wildcards.collapsed),
            'index': 'bowtie/{index}.1.ebwt'.format(index=index_name)}

# map all the reads to proper bowtie index
rule map_reads:
    input:
        unpack(map_reads_helper)
    output:
        sam=temp("{sample}/tmp/{sample}_{collapsed}.extracted.trim.sam"),
        log="{sample}/log/{sample}_{collapsed}.bowtie.log"
    params:
        index=lambda wildcards: 'bowtie/{}'.format(sample_sheet.loc[sample_sheet.ID == sample_name_to_id(wildcards.sample), 'ProbeName'].values[0])
    shell:
        "bowtie {params.index} <( gunzip < {input.trim} ) -v {mismatches} --sam > {output.sam} 2> {output.log}"

# creates a bowtie index for each sample (if there are different ones)
rule create_bowtie_index:
    input:
        lambda wildcards: sample_sheet.loc[sample_sheet.ProbeName == wildcards.index, 'ProbePath'].values[0]
    output:
        'bowtie/{index}.1.ebwt'
    params:
        import_mips=path.join(mip_directory, 'import_mips.py'),
        index=lambda wildcards: 'bowtie/{index}'.format(index=wildcards.index)
    shell:
        '{params.import_mips} -i {input} -o {params.index}.fa ; bowtie-build {params.index}.fa {params.index}'

# trim remaining read to get just 25bp homology region
# NOTE: this trimming isn't specific for bulk/emul, but we need to tell extract which one to do
rule trim_read:
    input:
        lambda wildcards: "{sample}/tmp/{sample}_R1_001.bulk_extracted.fastq.gz".format(sample=wildcards.sample) if sample_sheet.loc[sample_sheet.ID == sample_name_to_id(wildcards.sample), 'Read2'].values[0] == 'Bulk' else "{sample}/tmp/{sample}_{collapsed}_R1_001.emul_extracted.fastq.gz".format(sample=wildcards.sample, collapsed=wildcards.collapsed)
    output: 
        temp("{sample}/tmp/{sample}_{collapsed}_R1_001.extracted.trim.fastq.gz")
    shell:
        "zcat {input} | fastx_trimmer -o {output} -f 1 -l {readlen} -z -Q33"

# extraction rules (experiment-type specific)

rule get_bulk_extracted_reads:
    input:
        "{sample}/tmp/{sample}_R1_001.bulk_extracted.fastq.gz"
    output:
        "{sample}/stats/{sample}_bulk_extracted.stats.txt"
    shell:
        "zcat {input} | echo -e \"ExtractedReads\t$((`wc -l`/4))\" > {output}"

# extract umis from reads
rule bulk_extract:
    input:
        # bulk_extract_helper
        "{sample}/tmp/{sample}_R1_001.downsampled.fastq.gz"
    output:
        r1=temp("{sample}/tmp/{sample}_R1_001.bulk_extracted.fastq.gz"),
        log="{sample}/log/{sample}.extract.log"
    shell:
        "umi_tools extract --stdin {input} --stdout {output.r1}  \
            --bc-pattern={umi} -L {output.log}"

rule get_emul_extracted_reads:
    input:
        "{sample}/tmp/{sample}_{collapsed}_R1_001.emul_extracted.fastq.gz"
    output:
        "{sample}/stats/{sample}_emul_extracted_{collapsed}.stats.txt"
    shell:
        "zcat {input} | echo -e \"ExtractedReads_{wildcards.collapsed}\t$((`wc -l`/4))\" > {output}"

# extract umis/cbcs from reads (that match valid cbcs)
rule emul_extraction:
    input:
        r1='{sample}/tmp/{sample}_R1_001.filtered.fastq.gz',
        r2='{sample}/tmp/{sample}_R2_001.filtered.fastq.gz',
        wl='{sample}/tmp/{sample}_{collapsed}_whitelist.txt'
    output:
        r1=temp("{sample}/tmp/{sample}_{collapsed}_R1_001.emul_extracted.fastq.gz"),
        log="{sample}/log/{sample}_{collapsed}.extract.log"
    params:
        error_correct=lambda wildcards: '--error-correct-cell' if '_collapsed' in wildcards.collapsed else ''
    shell:
        "umi_tools extract --stdin {input.r1} --read2-in={input.r2} \
            --stdout {output.r1} --bc-pattern={umi} --bc-pattern2={cbc} \
            -L {output.log} --filter-cell-barcode --whitelist={input.wl} {params.error_correct}" #" --error-correct-cell" 

# grab starting cell barcodes (before collapsing) to write to file
rule get_starting_cells:
    input:
        "{sample}/tmp/{sample}_uncollapsed_whitelist.txt" # just to wait until whitelisting is over
    output:
        starting_cells="{sample}/stats/{sample}_starting_cells.stats.txt"
    params:
        guess=lambda wildcards: sample_sheet.loc[sample_sheet.ID == sample_name_to_id(wildcards.sample), 'CellNum'].values[0]
    run:
        guess = params.guess
        with open(output.starting_cells, 'w') as output:
            if guess == 'Guess':
                with open('{}/{}_cell_thresholds.tsv'.format(wildcards.sample, wildcards.sample)) as thresholds:
                    starting_bcs = 1 + int(thresholds.readlines()[1].strip())
            else:
                starting_bcs = guess
            output.write('StartingBCs\t{}\n'.format(starting_bcs))

# figures out how to call UMITools to determine cell number
def cell_number_from_samplesheet(wildcards):
    sample = wildcards.sample
    sample_id = sample_name_to_id(sample)
    cell_num = sample_sheet.loc[sample_sheet.ID == sample_id, 'CellNum'].values[0]
    if cell_num == 'Guess':
        return "--knee-method=distance"
    else:
        return "--set-cell-number={}".format(int(cell_num))

# select only valid cell barcodes
rule whitelist_cbcs:
    input:
        r1='{sample}/tmp/{sample}_R1_001.filtered.fastq.gz',
        r2='{sample}/tmp/{sample}_R2_001.filtered.fastq.gz'
    output:
        wl="{sample}/tmp/{sample}_uncollapsed_whitelist.txt",
        log="{sample}/log/{sample}.whitelist.log"
    params:
        plot_params="{sample}/{sample}",
        cell_num=cell_number_from_samplesheet
    shell:
        "umi_tools whitelist --stdin {input.r1} --read2-in={input.r2} \
            --bc-pattern={umi} --bc-pattern2={cbc} \
            --plot-prefix={params.plot_params} --method=umis \
            -L {output.log} {params.cell_num} --error-correct-threshold 0 > {output.wl}"

rule get_filtered_reads:
    input:
        "{sample}/tmp/{sample}_R1_001.filtered.fastq.gz"
    output:
        "{sample}/stats/{sample}_filtered.stats.txt"
    shell:
        "zcat {input} | echo -e \"FilteredReads\t$((`wc -l`/4))\" > {output}"

rule filter_fastqs:
    input:
        r1="{sample}/tmp/{sample}_R1_001.downsampled.fastq.gz",
        r2="{sample}/tmp/{sample}_R2_001.downsampled.fastq.gz"
    output:
        r1=temp("{sample}/tmp/{sample}_R1_001.filtered.fastq.gz"),
        r2=temp("{sample}/tmp/{sample}_R2_001.filtered.fastq.gz"),
        html="{sample}/log/{sample}_fastp.html",
        json="{sample}/log/{sample}_fastp.json"
    params:
        filter=path.join(mip_directory, 'filter_paired_fastqs_by_length.py'),
        fastp='/seq/lincRNA/Projects/MIP/jesse/190509_Test/fastp'
    shell:
        "{params.fastp} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            --length_required {cbc_len} -h {output.html} -j {output.json}"

rule get_downsampled_reads:
    input:
        "{sample}/tmp/{sample}_R1_001.downsampled.fastq.gz"
    output:
        "{sample}/stats/{sample}_downsampled.stats.txt"
    shell:
        "zcat {input} | echo -e \"InputReads\t$((`wc -l`/4))\" > {output}"

def downsample_helper(wildcards):
    sample, read = wildcards.sample, wildcards.read
    fastq = sample_sheet.loc[sample_sheet.ID == sample_name_to_id(sample), 'Read{}'.format(read)].values[0]
    return fastq

def get_read_count_for_downsampling(wildcards):
    read_count = sample_sheet.loc[sample_sheet.ID == sample_name_to_id(wildcards.sample), 'ReadsToUse'].values[0]
    return int(4 * read_count)

# downsample a fastq (for bulk or emulsion) by taking the first n reads (from sample sheet)
# returns all the reads if ReadsToUse is False
rule downsample:
    input:
        downsample_helper
    output:
        temp("{sample}/tmp/{sample}_R{read}_001.downsampled.fastq.gz")
    params:
        read_count=get_read_count_for_downsampling
    run:
        if params.read_count:
            # take the top n reads
            shell("set +o pipefail ; zcat {input} | head -n {params.read_count} | gzip > {output}")
        else:
            # take everything
            shell("cp {input} {output}")
            



