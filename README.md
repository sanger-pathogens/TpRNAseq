# Tp RNAseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**Tp RNAseq** is a pipeline for quantifying RNAseq reads that map to features of a bacterial genome using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [HTSeq](https://htseq.readthedocs.io). It has been adapted for use with *Treponema pallidum*, but could also in principle be used to analyse data for other bacterial species.

The pipeline can be run on local machines and high-performance computing (HPC) clusters with minimal setup.

## Pipeline summary

![TpRNAseq_workflow](./docs/images/TpRNAseq_workflow.svg)

In its simplest usage, **Tp RNAseq** takes a sample manifest (CSV; see [Generating a manifest](#generating-a-manifest)), reference (fasta) and annotation (gff) as input. It first runs some basic QC using `fastp` on the `.fastq.gz` files provided in the sample manifest. This involves adaptor removal and trimming of poor quality bases from the ends of the reads. The pipeline will then map reads to the given reference genome using Bowtie2. In order to do this, bowtie2 index files are required. These will be created if necessary.

Following mapping, the pipeline will mark and remove duplicates using Picard's [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard). These deduplicated alignments will then be quantified per feature using HTSeq. For each sample, the pipeline outputs fastqc reports pre- and post-QC and a HTSeq count table. The count tables are combined across samples to generate a summary count table for convenient downstream analysis. It also allows the user to optionally output alignment files (bam) pre- and post-deduplication. Reports are summarized using [multiQC](https://multiqc.info/).

## Requirements

### Software

- A Unix-like operating system environment (e.g. Linux, macOS, Windows with [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)) with Bash 3.2 or later.
- Java 11 or later ([OpenJDK](https://openjdk.org/)) as a nextflow dependency
- [nextflow](https://www.nextflow.io/) as workflow management system
- [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/), or [conda]() for pipeline dependency management

### Hardware (recommended)

- \>= 16GB RAM
- \>= 2 CPUs/cores
- Enough space for images (< 10GB) and intermediate data

## Dependency Installation
- **Linux / macOS**
    - Install [Java and Nextflow](https://www.nextflow.io/docs/latest/install.html)
    - Install any of:
        - [Docker Desktop](https://docs.docker.com/desktop/install) (recommended)
        - [Docker Engine](https://docs.docker.com/engine/install)
        - [Singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html)
        - [Apptainer](https://apptainer.org/)
        - [miniconda3](https://docs.anaconda.com/free/miniconda/)

- **Windows (with WSL)**
    - Install [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) by Microsoft
    - Install [Java and Nextflow](https://www.nextflow.io/docs/latest/install.html)
    - Install any of:
        - [Docker Desktop with WSL integration](https://docs.docker.com/desktop/wsl/)
        - [miniconda3](https://docs.anaconda.com/free/miniconda/)

- Troubleshooting:
    - After installation of Java, you may need to `export JAVA_ARG=$(which java)` to ensure successful nextflow installation.
    - Docker will usually require administrator privileges, so you may have to seek approval if using a managed machine (better to use `conda` in this case).
    - After installation of Docker Desktop on macOS, you might need to [allow Docker to access more system resources](https://docs.docker.com/desktop/settings/mac/), especially CPU and Memory, to match pipeline hardware requirements.

## Getting started

### Running on a personal computer

1.  [Install the necessary dependencies](#dependency-installation)

2.  Open a terminal.

3.  Run the pipeline:  
    For example input, please see [Generating a manifest](#generating-a-manifest).
    
    To run the pipeline with `docker`, use the `-profile docker` option:
    ```
    nextflow run sanger-pathogens/TpRNAseq \
        -profile docker \
        --manifest test_data/manifest.csv \
        --reference test_data/ref.fasta \
        --annotation test_data/ref.gff \
        --library_strandedness reverse
    ```

    Other profiles are also supported (`docker`, `conda`, `singularity`, `apptainer`).  
    :warning: If no profile is specified the pipeline will run with a Sanger HPC-specific configuration.

    See [usage](#usage) for all available pipeline options.

### Running on the farm (Sanger HPC clusters)

1.  Load nextflow and singularity modules:
    ```bash
    module load nextflow ISG/singularity
    ```

2.  Clone the repo and change directory into it.

3.  Start the pipeline:  
    For example input, please see [Generating a manifest](#generating-a-manifest).

    Example:
    ```bash
    nextflow run main.nf --manifest test_data/manifest.csv --reference test_data/ref.fasta --annotation test_data/ref.gff --library_strandedness reverse
    ```

    It is good practice to submit a dedicated job for the nextflow master process (use the `oversubscribed` queue):
    ```bash
    bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000 \
        nextflow run main.nf --manifest test_data/manifest.csv --reference test_data/ref.fasta --annotation test_data/ref.gff --library_strandedness reverse
    ```

    See [usage](#usage) for all available pipeline options.

4.  Once your run has finished, check output in the directory supplied to the `--outdir` option (default `./results`). Remember to clean up any intermediate files. To do this (assuming no other pipelines are running from the current working directory) run:

    ```bash
    rm -rf work .nextflow*
    ```

    Alternatively, use [`nextflow clean`](https://www.nextflow.io/docs/latest/cli.html#clean), which allows selective removal of these files.

### Running on any other HPC cluster

This pipeline should be compatible with most HPC clusters, but will not have the necessary configuration files for those clusters hosted in this repository.

Nevertheless, your insitution may support an [nf-core configuration profile](https://github.com/nf-core/configs/tree/master/conf). Documentation for the use of these profiles can be found [here](https://github.com/nf-core/configs?tab=readme-ov-file#documentation).

These profiles are supported out-of-the-box by this pipeline. You simply have to run the pipeline with the relevant profile. For example, to run the pipeline on the Cambridge HPC, run with `-profile cambridge`.

:warning: If no profile is specified the pipeline will run with a Sanger HPC-specific configuration. 


## Generating a manifest

Manifests supplied as an argument to `--manifest`, should be of of the following format:

```console
ID,REP,R1,R2
sample,1,./test_data/inputs/sample_rep1_1.fastq.gz,./test_data/inputs/sample_rep2_2.fastq.gz
sample,2,./test_data/inputs/sample_rep2_1.fastq.gz,./test_data/inputs/sample_rep2_2.fastq.gz
```

Where column `ID` can be an arbitrary sample identifier, `REP` describes the replicate structure of the data, `R1` is a .fastq.gz file of forward reads, `R2` is the mate .fastq.gz file containing reverse reads. The above example shows entries for two replicates of one sample.

## Usage

```console
---------------------------------------------------------------------------------------


88888888888            8888888b.  888b    888        d8888                           
    888                888   Y88b 8888b   888       d88888                           
    888                888    888 88888b  888      d88P888                           
    888  88888b.       888   d88P 888Y88b 888     d88P 888 .d8888b   .d88b.   .d88888
    888  888 "88b      8888888P"  888 Y88b888    d88P  888 88K      d8P  Y8b d88" 888
    888  888  888      888 T88b   888  Y88888   d88P   888 "Y8888b. 88888888 888  888
    888  888 d88P      888  T88b  888   Y8888  d8888888888      X88 Y8b.     Y88b 888
    888  88888P"       888   T88b 888    Y888 d88P     888  88888P'  "Y8888   "Y88888
         888                                                                      888
         888                                                                      888
         888                                                                      888
                              ____          ,--r ~.       ,^"   `"w                   
               ____        x^      'w     |L.._    ^m    A^''w     V                  
              D    "W     [R``'w    '@   ,R   [.    [L  jR    K     K     m" ` `W     
             A0     [H    R    [H    [   R     0     @  R     [     [L   0      R     
            ["[L     0   R     [@     @ R      RH    [ A       K     W  A      A      
           #"  [      K,R     /L[     0R     ,R D    [R       #0     [ #      #       
         y^    ,%     '@    ,R   L    [     z"  [    [H      A [L     R      A       
          `  "   T_    'W ,^     0    [L_,.#     L   [H    zC   0     @     R         
                  ^~-.. <^        %_     ,4`     T,   " ^"x"    !L    0__.gR         
                                    ``" `          "----^"        Y.____x^"           
---------------------------------------------------------------------------------------

 Usage: 

      nextflow run main.nf --manifest <manifest> --annotation <gff> --reference <fasta> --library_strandedness [reverse] --outdir [./results]

---------------------------------------------------------------------------------------
 Input 
      --manifest
            Path to a CSV manifest comprising 4 columns (ID, REP, R1, R2). ID is an arbitrary sample ID, REP describes replicate structure, R1 and R2 columns contain paths to *.fastq.gz files.
            Default: none
      --reference
            Path to reference genome in FASTA format.
            Default: none
      --annotation
            Path to genome annotation in GFF format.
            Default: none
      --library_strandedness
            Strandedness of the RNAseq library. Options: reverse, forward, none.
            Default: reverse
---------------------------------------------------------------------------------------

 Output 
      --outdir
            Path to an output directory (created if it doesn't exist).
            Default: ./results
      --keep_combined_fastqs
            If true, publish the combined fastqs.
            Default: false
      --keep_sorted_bam
            If true, keep the sorted bam file generated from mapping.
            Default: false
      --keep_dedup_bam
            If true, keep the deduplicated bam file.
            Default: false
      --keep_filtered_bam
            If true, keep filtered bam file.
            Default: false
---------------------------------------------------------------------------------------

 Data Combining 
      --combine_fastqs
            If true, combine fastqs at replicate level, based on sample identifiers provided in the manifest.
            Default: false
      --combine_rep
            If true, combine fastqs at the sample level, based on the replicate structure described in the manifest.
            Default: false
---------------------------------------------------------------------------------------

 QC 
      --trimmer
            Software to use for trimming. Options: fastp.
            Default: fastp
      --fastp_args
            Options and arguments that will be supplied to fastp to modify QC behaviour.
            Default: none
      --multiqc_config
            A multiqc configuration file (yml) that can be used to customise the multiqc behaviour.
            See https://multiqc.info/docs/getting_started/config/.
            Default: none
---------------------------------------------------------------------------------------

 Mapping 
      --bowtie2_args
            Options and arguments that will be supplied to Bowtie2 to modify mapping behaviour.
            Default: --local --very-sensitive-local --rdg 8,4 --rfg 8,4 --no-mixed
      --dedup
            If true, remove duplicates (hopefully library prep artifacts) from bam file.
            Default: false
---------------------------------------------------------------------------------------

 Counting 
      --htseq_args
            Options and arguments that will be supplied to htseq-count to modify counting behaviour.
            Default: --type gene --idattr locus_tag --nonunique none --secondary-alignments ignore
      --samtools_filter_args
            Arguments supplied to samtools to will be used to filter alignments of interest for counting. Acceptable arguments: -f or -F.
            The default will only keep reads that aligned in proper pairs.
            Default: -f 2
---------------------------------------------------------------------------------------

 Coverage 
      --strand_specific
            If true, run strand-specific coverage analysis.
            Default: false
      --coverage_window_size
            Size of window or step over which to compute cumulative per base coverage (output as a wig file).
            Default: 100
      --coverage_context
            Provide context around the annotated region (coverage_window_size base pairs either side) in the coverage plot.
            Default: 100
      --pairwise
            If true, generate pairwise sample comparison coverage plots.
            Default: false
---------------------------------------------------------------------------------------

 Logging 
      --monochrome_logs
            Should logs appear in plain (non-coloured) ASCII
            Default: false
---------------------------------------------------------------------------------------
```

## Output

The pipeline will generate a directory structure in the output directory specified by `--outdir`. Some directories may be missing depending on the data retention flags employed (See [Usage](#usage)). Note that parameter values are represented by `{param_name}` below.
```
results/
├── bowtie2
├── coverage
│   ├── plots
│   ├── wig_{coverage_window_size}
│   └── wig_raw
├── fastp
│   └── fastqc
│       ├── raw
│       └── trim
├── filtered_bams
│   ├── minus_filter
│   ├── plus_filter
│   └── user_defined_filter
├── htseq
├── multiqc
│   ├── multiqc_report_data
│   └── multiqc_report.html
├── pipeline_info
├── sorted_bams
└── sorted_ref
```

The content of these folders is described in more detail below. Note that output files will be sample-specific if `--combine_rep` is used and replicate-specific otherwise:
| Directory / File | Description |
| --- | ---|
| `bowtie2` | Contains Bowtie2 index files (`.bt2.*`) for the given reference (`--reference`). |
| `coverage/plots` | Contains per sample (or pairwise-sample if `--pairwise`) folders, containing strand-specific coverage plots for all `gene` level annotations in the given `--gff` and intergenic regions. |
| `coverage/wig_{coverage_window_size}` | Contains `*.wig` files which summarize cumulative (per base) genome coverage over steps of the given `--coverage_window_size`. |
| `coverage/wig_raw` | Contains `*.wig` files per base genome coverage. |
| `fastp/trimmed_fastqs` | Contains trimmed reads in `*.fastq.gz` files.<br> Generated with `--keep_trimmed_fastqs`. |
| `fastp/reports` | Contains fastp report files in `html` and `json` format. |
| `fastp/fastqc` | Contains FASTQC data and reports for both trimmed and raw `*.fastq.gz` files. |
| `filtered_bams/(plus\|minus)_filter` | Contains reverse (`plus_filter`) and forward (`minus_filter`) strand-specific BAM files, containing reads that mapped to each strand of the given `--reference`.<br> Generated with `--keep_filter_bam`. |
| `filtered_bams/user_defined_filter` | Contains BAM files, containing reads that mapped to the given `--reference` and filtered using `--samtools_filter_args`.<br> Generated with `--keep_filter_bam`. |
| `htseq` | Contains HTSeq-generated count data (`tsv` format) suitable for differential expression analysis. It also contains a summary across all samples `gene_counts.tsv`. |
| `multiqc` | Contains a summary of all interesting QC metrics/statistics generated by the pipeline. |
| `sorted_bams` | Contains raw sorted BAM files (prior to filtering).<br> Generated with `--keep_sorted_bam`. |
| `sorted_ref` | Contains a copy of the reference `fasta` and fasta index `*.fai`.<br> Generated only if a reference `*.fai` could not be detected in the same directory as the input reference. |
| `pipeline_info` | Contains an execution report (`html`), execution trace (`tsv`) - useful for identifying failed processes, and a pipeline flowchart (`svg`). |


## Credits

**Tp RNAseq** was produced by PAM informatics at the Wellcome Sanger Institute, in collaboration with Linda Grillova and Eli Carrami. Inspiration was taken, in places, from the [nf-core pipelines](https://nf-co.re/pipelines).

## Support

If you require any help running this pipeline, experience a bug, or would like to request a new feature, please post an issue.

## Citations

This pipeline was developed to support the work of the following publications. If you use this pipeline in your work, please cite as appropriate:

<TODO - List of citations>

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
