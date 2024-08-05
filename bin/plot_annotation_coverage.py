#!/usr/bin/env python3

import argparse
import logging
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from matplotlib.ticker import FuncFormatter
from pathlib import Path
from typing import Generator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', '-s', dest='sample', help="Name of a sample for which to show coverage")
    parser.add_argument('--gff', '-g', dest="gff", help="Genome annotation in GFF3 format", type=Path)
    parser.add_argument('--wig_dir', '-d', dest="wig_dir", help="Path to directory containing .wig coverage files", type=Path)
    parser.add_argument('--ext', '-e', dest="ext", help="Number of bases by which to extend coverage plot either side of an annotated gene", type=int, default=100)
    parser.add_argument('--min_intergenic', '-i', dest="min_intergenic", help="Minimum number of bases for an intergenic region to be annotated", type=int, default=10)
    parser.add_argument('--outdir', '-o', dest="outdir", help="Path to an output directory in which to save the coverage plots", type=Path, default="./plots")
    return parser.parse_args()


def get_gene_name(ann_row: pd.DataFrame):
    """
    Extract gene name and locus identifier from attributes (tags) column origingating from GFF

    For locus id, prefer old_locus_tag to locus_tag attr
    For gene name, prefer gene attr, use locus id otherwise
    """
    t = ann_row["tags"].split(";")
    d = {i.split("=")[0]: i.split("=")[1] for i in t}
    if "old_locus_tag" in d:
        lt = d["old_locus_tag"]
    else:
        lt = d["locus_tag"]
    if "gene" in d:
        return d["gene"], lt
    else:
        return lt, lt


def load_data(sample_id: str, data_dir: Path, normalize: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load wig files, extract the coverage column, normalise the coverage, collect into dfs for plus and minus strand"""
    plus_data, minus_data = [], []
    for i, f in enumerate(data_dir.glob(f"{sample_id}_REP*plus.wig"), start=1):
        # Extract coverage (and label replicate R1, R2, etc.)
        p_data = pd.read_csv(f, sep="\t", header=None, names=['chrom', 'pos', 'cov'])['cov']
        #TODO fix potential bug if replacing plus with minus in full path (rather than just filename)!
        m_data = pd.read_csv(str(f).replace('plus', 'minus'), sep="\t", header=None, names=['chrom', 'pos', 'cov'])['cov']
        p_data.name = f'rep{i}'
        m_data.name = f'rep{i}'
        # Normalise
        if normalize:
            #TODO Check if this normalization is appropriate!
            nf = p_data.sum() + m_data.sum()
            plus_data.append(1e6 * p_data / nf)
            minus_data.append(1e6 * m_data / nf)
        else:
            plus_data.append(p_data)
            minus_data.append(m_data)
    return pd.concat(plus_data, axis=1), pd.concat(minus_data, axis=1)


def get_gene_annotations(gff: Path) -> pd.DataFrame:
    ann_base = pd.read_csv(gff, sep="\t", skiprows=3, header=None)
    # Get gene entries from GFF
    ann_base = ann_base[ann_base[2] == 'gene']
    # Subset "start", "end", "strand", "tags"
    ann_base = ann_base[[3, 4, 6, 8]]
    ann_base.columns = ["start", "end", "strand", "tags"]
    return ann_base

def add_intergenic_region_annotations(ann_base: pd.DataFrame, min_length: int = 10) -> list:
    IG_idx = 1
    ann = []
    # num_genes_in_gff = range(len(ann_base))
    for idx in range(len(ann_base)):
        name = get_gene_name(ann_base.iloc[idx])
        prev_name = get_gene_name(ann_base.iloc[idx-1])
        # If previous "gene" entry in GFF ends more than min_length bases before the current entry start
        # Determines the minimum length of integergenic region to be annotated
        # TODO Has the first IG region been accounted for in this code? Not really, in the sense that if a IG region exists between last and first gene, and (start of the first gene - min_length) < 0 then no matter what the end of the previous gene is, the if statement will never be satisfied and the IG region will never be added. BUT... it is added in an ad hoc manner at the end of this function!
        # Assume no overlapping genes with start of first gene. Note that normally the start of coordinates should fall in the intergenic region before the replication site (origin or replication - ori? or DnaA gene?).
        if ann_base.iloc[idx-1]["end"] < ann_base.iloc[idx]["start"] - min_length:
            # Then add an intergenic annotation
            ann.append([
                f"IG_{IG_idx}",  # locus id
                f"Inter_{prev_name[0]}:{name[0]}",  # gene name (or more human-meaningful name)
                ann_base.iloc[idx-1]["end"] + 1,  # start
                ann_base.iloc[idx]["start"] - 1,  # end
                '0'  # strand
            ])
            IG_idx += 1
        # Append the current entry annotation
        ann.append([name[1], name[0], ann_base.iloc[idx]["start"], ann_base.iloc[idx]["end"], ann_base.iloc[idx]["strand"]])
    # Append the last IG region (assumes there is one)
    ann.append([
        f"IG_{IG_idx}",
        f"Inter_{name[0]}:{get_gene_name(ann_base.iloc[0])[0]}",
        ann_base.iloc[idx]["end"] + 1,
        ann_base.iloc[0]["start"] - 1,
        '0'
    ])
    return ann
    

def generate_gene_with_neighbours(ann: list) -> Generator:
    for i in range(len(ann)):
        ## Identifying the previous (h) and next (j) gene
        # Special case: Last annotation
        if i == len(ann) - 1:
            # Since T pallidum genome is circular, next gene is first gene
            h, j = i-1, 0
        else:
            # Normal case
            h, j = i-1, i+1
        # If prev "gene" would be an intergenic region:
        if ann[h][4] == '0':
            # Subtract again, to get true previous gene
            h -= 1
        # If next "gene" would be an intergenic region:
        if ann[j][4] == '0':
            # Add again, to get true next gene
            j += 1
        # If the next "gene" would go out of the range, i.e. equal to len(ann), next gene = first gene
        # I think this should be handled by the special case above.
        if j == len(ann):
            j = 0
        yield ann[h], ann[i], ann[j]  # prev, gene, next


def get_coverage_context(coverage_data: pd.DataFrame, gene, chr_len):
    if gene[3] < gene[2]:
        gene[3] += chr_len
    start, end = get_region_limits(gene)
    if start < 0:
        new_cov_data = np.concatenate((coverage_data[start:], coverage_data[:end]))
    elif end > chr_len:
        new_cov_data = np.concatenate((coverage_data[start:], coverage_data[:end-chr_len]))
    else:
        new_cov_data = coverage_data[start:end]
    return new_cov_data


def get_region_limits(gene: list, ext: int = 100) -> tuple[int, int]:
    """
    Get the plotting limits for a regions corresponding to the given gene

    Args:
        gene (list): List of properties for the gene/intergenic region (locus id, name, start, end, strand)
        chr_len (int): Length of chromosome
        ext (int, optional): Region either side of gene to consider in plot, helpful due to uncertainty in 5' and 3' UTR. Defaults to 100.

    Returns:
        tuple[int, int]: Limits of the region (for plotting)
    """
    #TODO Original code had start = (gene[2] - ext - 1). Was that better?
    return (gene[2] - ext), (gene[3] + ext)


def get_chr_len(raw_coverage: pd.DataFrame) -> int:
    return len(raw_coverage)


def get_fragment_size(gene: list) -> int:
    """Size of a gene or feature"""
    return gene[3] - gene[2] + 1  # + 1 because GFF coordinates are inclusive


def plot_coverage(axis: matplotlib.axes.Axes, coverage_data: pd.DataFrame, gene: list, line_style: str = '-') -> matplotlib.axes.Axes:
    start, end = get_region_limits(gene)
    axis.set_xlim(start, end)
    axis.plot(range(start, end), coverage_data, line_style)
    return axis


def mark_chr_len(gene: list, chr_len: int):
    start, end = get_region_limits(gene)
    if end > chr_len:
        plt.axvline(x=chr_len, color='k', lw=2, ls="-")


def mark_gene_boundaries(gene: list) -> None:
    """
    Mutates currently active plot with vertical lines denoting annotated start and end of gene
    """
    gene_start, gene_end = gene[2], gene[3]
    plt.axvline(x=gene_start, color='grey', lw=1, ls="--")
    plt.axvline(x=gene_end, color='grey', lw=1, ls="--")
    

def highlight_neighbouring_genes(axis: matplotlib.axes.Axes, gene_with_neighbours: tuple[list, list, list], chr_len: int,  colours: tuple[str, str] = ('b','r')):
    prev, gene, next = gene_with_neighbours
    start, end = get_region_limits(gene)
    max_val = axis.get_ylim()[1]
    sense_col, anti_col = colours
    left_col = anti_col if prev[4] == '-' else sense_col
    right_col = anti_col if next[4] == '-' else sense_col
    prev_end = prev[3] - chr_len if start < 0 else prev[3]
    next_start = next[2] + chr_len if next[2] < gene[3] else next[2]
    if prev_end > start:
        axis.hlines(y=max_val, xmin=start, xmax=prev_end + 1, color=left_col, lw=10, alpha=0.5)
    if next_start < end:
        axis.hlines(y=max_val, xmin=next_start, xmax=end, color=right_col, lw=10, alpha=0.5)
    return axis


def add_gene_arrows(ax1: matplotlib.axes.Axes, gene: list, colours: tuple[str, str] = ('b', 'r')):
    # Because the max value changes if we add something near the top of the plot, the ylim has changed from what it was previously.
    # This causes next/prev gene annotations (if present) to appear below the (current) gene's arrow, but the plot still looks good!
    sense_col, anti_col = colours
    max_val = ax1.get_ylim()[1]
    start, end = get_region_limits(gene)
    region_size = end - start
    fragment_size = get_fragment_size(gene)
    arrow_kwargs = {
            "width": max_val / 50,  # simply 50x less than the height of the plot
            "length_includes_head": True,
            "head_length": region_size / 25,
        }
    if gene[4] == '+':
        ax1.arrow(gene[2], max_val, fragment_size, 0,  facecolor=sense_col, **arrow_kwargs)
    elif gene[4] == '-':
        ax1.arrow(gene[3], max_val, -fragment_size, 0, facecolor=anti_col, **arrow_kwargs)

def add_plot_title(gene):
    fragment_size = get_fragment_size(gene)
    strand = {"+": "+ strand", "-": "- strand", "0": ""}
        # If there's no readable gene name
    if gene[1] == gene[0]:
        plt.title(f"{gene[1]} [{fragment_size} bp] {strand[gene[4]]}")
    else:
        # If there's a readable gene name
        plt.title(f"{gene[1]} ({gene[0]}) [{fragment_size} bp] {strand[gene[4]]}")

def update_x_axis(axis: matplotlib.axes.Axes, chr_len: int):
    start, end = axis.get_xlim()
    region_size = end - start
    print(region_size)
    if region_size < 400:
        step_size = 100 * (region_size // 200)
    else:
        step_size = 100 * (region_size // 400)
    # start from the start of the region rounded up to the nearest hundred, cover the whole region, in step sizes described above
    ticks = np.arange(100 * (start // 100 + 1), end, step_size)
    axis.xaxis.set_ticks(ticks)
    if start >= 10**6:
        # Update font size (or could rotate) to avoid label clash
        axis.set_xticklabels(ticks, fontsize=8)
        # axis.set_xticklabels(ticks, rotation=45, ha='right')
    axis.xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x % chr_len), ',')))


def write_plot_to_file(save_path: Path, gene: list, i: int, num_ann: int):
    zfill_size = int(np.ceil(np.log10(num_ann)))
    readable_gene_name = gene[1]
    plt.savefig(f"{save_path}/feature_{str(i).zfill(zfill_size)}_{readable_gene_name}.png")


def main():
    args = parse_args()

    save_path = args.outdir / f"{args.sample}"  #TODO Keep sample pair named subdirectory if we now give the user the choice of where to save? Probably not...
    os.makedirs(save_path, exist_ok=True)

    # Parse wig files for given sample into dfs per base coverage data
    p1, m1 = load_data(args.sample, args.wig_dir)
    chr_len = get_chr_len(p1)

    # Parse GFF
    ann_base = get_gene_annotations(args.gff)
    ann = add_intergenic_region_annotations(ann_base, args.min_intergenic)

    # Loop over genes and plot
    gene_with_neighbours = list(generate_gene_with_neighbours(ann))
    for i, (prev, gene, next) in enumerate(gene_with_neighbours, start=1):
        logging.info(f"Processing {gene}") 
        sense_cov = get_coverage_context(p1, gene, chr_len)
        anti_cov = get_coverage_context(m1, gene, chr_len)

        # plot
        fig, ax1 = plt.subplots()

        plot_coverage(ax1, sense_cov, gene, line_style='b-')
        plot_coverage(ax1, anti_cov, gene, line_style='r-')
        mark_gene_boundaries(gene)
        mark_chr_len(gene, chr_len)
        highlight_neighbouring_genes(ax1, (prev, gene, next), chr_len)
        add_gene_arrows(ax1, gene)
        update_x_axis(ax1, chr_len)
        add_plot_title(gene)

        write_plot_to_file(save_path, gene, i, num_ann=len(ann))
        plt.clf()
        plt.close('all')

if __name__ == "__main__":
    main()
