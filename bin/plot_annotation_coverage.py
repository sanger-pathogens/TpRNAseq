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
    parser.add_argument('--sample_1', '-s1', dest='sample_1', help="Name of a sample for which to show coverage")
    parser.add_argument('--sample_2', '-s2', dest='sample_2', help="Name of a second sample for which to show coverage", required=False)
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


def load_data(name_id: str, data_dir: Path, normalize: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load wig files, extract the coverage column, normalise the coverage, collect into dfs for plus and minus strand"""
    plus_data, minus_data = [], []
    for i, f in enumerate(data_dir.glob(f"{name_id}*plus.wig"), start=1):
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

def add_coverage_data_context(coverage_data: dict[str, pd.DataFrame], ext: int):
    extended_coverage_data = {}
    for k, v in coverage_data.items():
        chr_len = len(v)
        v = v.reindex(range(1, chr_len + 1))
        print(v)
        lower_context = v[(chr_len - ext):(chr_len + 1)]
        lower_context = lower_context.reindex(range(-1 - ext, 0)) 
        upper_context = v[0:ext]
        upper_context.reindex(range(chr_len + 1, chr_len + ext + 1))
        extended_coverage_data[k] = pd.concat([lower_context, v, upper_context])
        print(extended_coverage_data[k][0:102])
        # print(extended_coverage_data[k][-101])
        # print(extended_coverage_data[k][1139616:1139717])


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

class feature:
    """Feature in genome, could be gene or intergenic region"""
    def __init__(self, id, name, start, end, strand, tags, prev, next):
        self.id
        self.name
        self.start
        self.end
        self.strand
        self.tags
        self.prev_gene
        self.next_gene


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
    return (gene[2] - ext), (gene[3] + ext)

def get_plot_data(coverage_data: dict[str, pd.DataFrame], region_limits: tuple) -> dict[str, dict]:
    """
    Generate data required for plotting coverage across loci in genome

    Args:
        coverage_data (dict): Dictionary mapping coverage data by strand_sample
        region_limits (tuple): Start and end of the region to be plotted
    """
    start, end = region_limits
    if len(coverage_data) == 2:
        plot_data = {
            'sense_1': {'color': 'b'}, 'anti_1': {'color': 'r'},
        }
    elif len(coverage_data) == 4:
        plot_data = {
            'sense_1': {'color': 'b'}, 'anti_1': {'color': 'r'},
            'sense_2': {'color': 'c'}, 'anti_2': {'color': 'm'},
        }
    else:
        raise ValueError(f"Unexpected number of elements in coverage data: ${coverage_data}")

    for k, v in coverage_data.items():
        # Find mean, min, max coverage across all replicates for this gene's region (on both strands, and in both samples)
        plot_data[k]['values'] = v.iloc[start:end].mean(1).to_numpy()
        plot_data[k]['min'] = v.iloc[start:end].min(1).to_numpy()
        plot_data[k]['max'] = v.iloc[start:end].max(1).to_numpy()
    return plot_data

def get_chr_len(raw_coverage: pd.DataFrame) -> int:
    return len(raw_coverage)

def get_fragment_size(gene: list) -> int:
    """Size of a gene or feature"""
    return gene[3] - gene[2] + 1  # + 1 because GFF coordinates are inclusive

def get_region_size(region_limits: tuple[int, int]) -> int:
    """Size of the region that is shown in a plot"""
    return region_limits[1] - region_limits[0]

def plot_base_line_graph(plot_data: dict, region_limits: tuple) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    fig, ax1 = plt.subplots()
    ax1.set_xlim(*region_limits)

    # Line graph of coverage
    for label, info in plot_data.items():
        x_range = range(*region_limits)
        ax1.plot(x_range, info['values'], info['color'] + '-')  # mean
        ax1.fill_between(x_range, info['min'], info['max'], color=info['color'], alpha=0.1)  # uncertainty
    return fig, ax1

def mark_gene_boundaries(gene: list) -> None:
    """
    Mutates currently active plot with vertical lines denoting annotated start and end of gene
    """
    plt.axvline(x=gene[2], color='grey', lw=1, ls="--")
    plt.axvline(x=gene[3], color='grey', lw=1, ls="--")

def highlight_neighbouring_genes(ax1: matplotlib.axes.Axes, gene_and_neighbours: tuple, chr_len: int):
    prev, gene, next = gene_and_neighbours
    # Assign overlap regions and their color
    # Choose colour for left and right ranges (100bp) either side of gene
    left_col = 'r' if prev[4] == '-' else 'b'
    right_col = 'r' if next[4] == '-' else 'b'
    max_val = ax1.get_ylim()[1]
    start, end = ax1.get_xlim()
    # left_box = coordinate of end of previous feature
    left_box = prev[3]
    # Draw horizontal lines at top of plot indicating whether next/prev feature is on + or - strand
    # Given that we extend the genes by 100bp, this handles the possibility of the extended region overlapping with the (unextended) range of the previous or next gene.
    if left_box > start:
        # if the the end of the prev gene is within the start of the (extended) region/plot
        # the highlight the from start of the plot to the end of the prev gene
        ax1.hlines(y=max_val, xmin=start, xmax=left_box+1, color=left_col, lw=10, alpha=0.5)
    if end > chr_len and next[2] < (end - chr_len):  # last gene
        # if end of the (extended) region/plot is beyond the end of the ref
        # and the start of the next gene (probably the first) is before the end of the extended region/plot (last gene),
        # then highlight between the start of the next gene and the end of the region/plot
        ax1.hlines(y=max_val, xmin=next[2]+chr_len, xmax=end, color=right_col, lw=10, alpha=0.5)  # next[2] should be 0?
        plt.axvline(x=chr_len, color='k', lw=2, ls="-")  # Mark the end of the ref with a solid line
    if next[2] < end and next[3] > gene[2]:
        # if the start of the next gene is within the end of the (extended) region/plot
        # and the end of the next gene is beyond the start of the current gene
        # then highlight between the start of the next gene and the end of the plot
        ax1.hlines(y=max_val, xmin=next[2], xmax=end, color=right_col, lw=10, alpha=0.5)

def add_gene_arrows(ax1: matplotlib.axes.Axes, gene: list, region_size: tuple[int], fragment_size: int):
    # Because the max value changes if we add something near the top of the plot, the ylim has changed from what it was previously.
    # This causes next/prev gene annotations (if present) to appear below the (current) gene's arrow, but the plot still looks good!
    max_val = ax1.get_ylim()[1]
    arrow_kwargs = {
            "width": max_val / 50,  # simply 50x less than the height of the plot
            "length_includes_head": True,
            "head_length": region_size / 25,
        }
    if gene[4] == '+':
        ax1.arrow(gene[2], max_val, fragment_size, 0,  facecolor='blue', **arrow_kwargs)
    elif gene[4] == '-':
        ax1.arrow(gene[3], max_val, -fragment_size, 0, facecolor='red', **arrow_kwargs)

def add_plot_title(gene, fragment_size):
    strand = {"+": "+ strand", "-": "- strand", "0": ""}
        # If there's no readable gene name
    if gene[1] == gene[0]:
        plt.title(f"{gene[1]} [{fragment_size} bp] {strand[gene[4]]}")
    else:
        # If there's a readable gene name
        plt.title(f"{gene[1]} ({gene[0]}) [{fragment_size} bp] {strand[gene[4]]}")

def update_x_axis(ax1: matplotlib.axes.Axes, chr_len: int, region_size: tuple[int, int]):
    start, end = ax1.get_xlim()
    if region_size < 400:
        step_size = 100 * (region_size // 200)
    else:
        step_size = 100 * (region_size // 400)
        # start from the start of the region rounded up to the nearest hundred, cover the whole region, in step sizes described above
    ax1.xaxis.set_ticks(np.arange(100 * (start // 100 + 1), end, step_size))
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x % chr_len), ',')))

def write_plot_to_file(save_path: Path, gene: list, i: int, num_ann: int):
    zfill_size = int(np.ceil(np.log10(num_ann)))
    readable_gene_name = gene[1]
    plt.savefig(f"{save_path}/feature_{str(i).zfill(zfill_size)}_{readable_gene_name}.png")

def contextualize_coordinates(start: int, end: int, chr_len: int):
    # Account for region where end overlaps with chromosome end
    # TODO This will never happen if we are simply extending the chromosome with context...
    if end < start:
        end += chr_len + 1  #TODO + 1 because 1-based indexing
    # Account for region where start overlaps with chromosome end
    if start < 0:
        start += chr_len + 1  #TODO + 1 because 1-based indexing
    return start, end

def contextualize_feature(feature: list, chr_len: int):
    feature = feature.copy()
    feature[2], feature[3] = contextualize_coordinates(feature[2], feature[3], chr_len)
    return feature

def contextualize_features(features: list, chr_len: int):
    new_features = []
    for feature in features:
        new_features.append(contextualize_feature(feature, chr_len))
    return new_features
    

def main():
    args = parse_args()

    # Parse wig files for given sample into dfs of normalised per base coverage count data
    if args.sample_2 is None:
        save_path = args.outdir / f"{args.sample_1}"  #TODO Keep sample pair named subdirectory if we now give the user the choice of where to save? Probably not...
        p1, m1 = load_data(args.sample_1, args.wig_dir)
        coverage_data = {
            'sense_1': p1,
            'anti_1': m1,
        }
    else:
        save_path = args.outdir / f"{args.sample_1}_{args.sample_2}"  #TODO Keep sample pair named subdirectory if we now give the user the choice of where to save? Probably not...
        p1, m1 = load_data(args.sample_1, args.wig_dir)
        p2, m2 = load_data(args.sample_2, args.wig_dir)
        coverage_data = {
            'sense_1': p1,
            'anti_1': m1,
            'sense_2': p2,
            'anti_2': m2
        }
    os.makedirs(save_path, exist_ok=True)
    chr_len = get_chr_len(p1)

    add_coverage_data_context(coverage_data, 100)

    # Parse GFF
    ann_base = get_gene_annotations(args.gff)
    ann = add_intergenic_region_annotations(ann_base, args.min_intergenic)
    # print(ann[-2], ann[-1], ann[0])

    # Loop over genes and plot
    gene_with_neighbours = generate_gene_with_neighbours(ann)
    gene_with_neighbours = list(gene_with_neighbours)
    print(gene_with_neighbours[0])
    # contextualize_coordinates(gene_with_neighbours, chr_len, args.ext)
    for i, (prev, gene, next) in enumerate(gene_with_neighbours, start=1):
        logging.info(f"Processing {gene}")
        prev, gene, next = contextualize_features([prev, gene, next], chr_len)
        print(prev, gene, next)
        region_limits = get_region_limits(gene, args.ext)
        # print(region_limits)
        region_limits = contextualize_coordinates(*region_limits, chr_len)
        print(region_limits)
        region_size = get_region_size(region_limits)
        # print(region_size)
        fragment_size = get_fragment_size(gene)
        # print(fragment_size)
        plot_data = get_plot_data(coverage_data, region_limits)

        # Plot comparison of coverage between the samples and strands
        fig, ax1 = plot_base_line_graph(plot_data, region_limits)
        mark_gene_boundaries(gene)
        highlight_neighbouring_genes(ax1, (prev, gene, next), chr_len)
        add_gene_arrows(ax1, gene, region_size, fragment_size)
        add_plot_title(gene, fragment_size)
        update_x_axis(ax1, chr_len, region_size)

        write_plot_to_file(save_path, gene, i, num_ann=len(ann))
        plt.clf()
        plt.close('all')
        break

if __name__ == "__main__":
    main()
