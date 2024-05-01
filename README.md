# Tp RNAseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**Tp RNAseq** is a pipeline for quantifying RNAseq reads that map to features of a bacterial genome using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [HTSeq](https://htseq.readthedocs.io). It has been adapted for use with *Treponema pallidum*, but could also in prinicple be used to analyse data for other bacterial species.

## Pipeline summary

In its simplest usage, **Tp RNAseq** takes a sample manifest (CSV; see [Generating a manifest](#generating-a-manifest)), reference (fasta) and annotation (gff) as input. It first runs some basic QC using `fastp` on the `.fastq.gz` files provided in the sample manifest. This involves adaptor removal and trimming of poor quality bases from the ends of the reads. The pipeline will then map reads to the given reference genome using Bowtie2. In order to do this, bowtie2 index files are required. These will be created if necessary.

Following mapping, the pipeline will mark and remove duplicates using Picard's [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard). These deduplicated alignments will then be quantified per feature using HTSeq. For each sample, the pipeline outputs fastqc reports pre- and post-QC and a HTSeq count table. The count tables are combined across samples to generate a summary count table for convenient downstream analysis. It also allows the user to optionally output alignment files (bam) pre- and post-deduplication. Reports are summarized using [multiQC](https://multiqc.info/).

## Getting started

### Running on the farm (Sanger HPC clusters)

1. Load nextflow and singularity modules:
   ```bash
   module load nextflow ISG/singularity
   ```

2. Clone the repo and change directory into it.

3. Start the pipeline:
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

4. Once your run has finished, check output in the directory supplied to the `--outdir` option (default `./results`). Remember to clean up any intermediate files. To do this (assuming no other pipelines are running from the current working directory) run:

   ```bash
   rm -rf work .nextflow*
   ```

   Alternatively, use [`nextflow clean`](https://www.nextflow.io/docs/latest/cli.html#clean), which allows selective removal of these files.

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
            Path to a CSV manifest comprising 3 columns (sample_id, R1, R2). R1 and R2 columns contain paths to *.fastq.gz files.
            Default: none
      --reference
            Path to reference genome in FASTA format.
            Default: none
      --annotation
            Path to genome annotation in GFF format.
            Default: none
      --library_strandedness
            Strandedness of the RNAseq library.
            Default: reverse
---------------------------------------------------------------------------------------

 QC 
      --trimmer
            Software to use for trimming. Options: fastp.
            Default: fastp
      --fastp_args
            Options and arguments that will be supplied to fastp to modify QC behaviour.
            Default: none
      --multiqc_config
            A multiqc configuration file (yml) that can be used to customise the multiqc behaviour. See https://multiqc.info/docs/getting_started/config/.
            Default: none
---------------------------------------------------------------------------------------

 Mapping 
      --bowtie2_args
            Options and arguments that will be supplied to Bowtie2 to modify mapping behaviour.
            Default: --local --very-sensitive-local --rdg 8,4 --rfg 8,4 --no-mixed
      --dedup
            If true, deduplicate bam file.
            Default: false
      --keep_sorted_bam
            If true, keep the sorted bam file generated from mapping.
            Default: false
      --keep_dedup_bam
            If true, keep the deduplicated bam file.
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

 Logging 
      --monochrome_logs
            Should logs appear in plain (non-coloured) ASCII
            Default: false
---------------------------------------------------------------------------------------
```

## Credits

**Tp RNAseq** was produced by PAM informatics at the Wellcome Sanger Institute.

## Support

If you require any help running this pipeline, experience a bug, or would like to request a new feature, please post an issue.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
