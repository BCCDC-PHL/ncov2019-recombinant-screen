# Recombinant Screen for SARS-CoV-2

This pipeline is based directly on [ktmeaton/ncov-recombinant](https://github.com/ktmeaton/ncov-recombinant).
We've adapted the original Snakemake pipeline into a Nextflow pipeine, and made minor adjustments to fit our
systems and SARS-CoV-2 sequence datasets.

## Overview

```mermaid
flowchart TD
  ref(ref.fa)
  primers(primer.bed)
  metadata(metadata.tsv)
  problematic_sites(problematic_sites.vcf)
  breakpoints(breakpoints.tsv)
  breakpoints(breakpoints.tsv) --> issues_download(issues_download)
  run_dir(run_dir) --> identify_complete_genomes(identify_complete_genomes)
  nextclade_dataset(nextclade_dataset)
  identify_complete_genomes(identify_complete_genomes) --> prepare_multi_fasta(prepare_multi_fasta)
  prepare_multi_fasta(prepare_multi_fasta) --> nextclade(nextclade)
  nextclade_dataset(nextclade_dataset) --> nextclade(nextclade)
  metadata(metadata.tsv) --> nextclade(nextclade)
  nextclade(nextclade) --> nextclade_recombinants(nextclade_recombinants)
  nextclade_recombinants(nextclade_recombinants) --> sc2rf(sc2rf)
  primers --> sc2rf(sc2rf)
  sc2rf(sc2rf) --> sc2rf_recombinants(sc2rf_recombinants)
  nextclade(nextclade) --> sc2rf(sc2rf)
  issues_download(issues_download) --> sc2rf(sc2rf)
  nextclade(nextclade) --> nextclade_exclude_false_positives(nextclade_exclude_false_positives)
  sc2rf_recombinants(sc2rf_recombinants) --> nextclade_exclude_false_positives(nextclade_exclude_false_positives)
  sc2rf_recombinants(sc2rf_recombinants) --> fasta_to_vcf(fasta_to_vcf)
  problematic_sites(problematic_sites.vcf) --> fasta_to_vcf(fasta_to_vcf)
  ref(ref.fa) --> fasta_to_vcf(fasta_to_vcf)
  usher_download(usher_download) --> usher(usher)
  usher_download(usher_download) --> usher_columns(usher_columns)
  fasta_to_vcf(fasta_to_vcf) --> usher(usher)
  usher(usher) --> usher_stats(usher_stats)
  issues_download(issues_download) --> usher_stats(usher_stats)
  metadata(metadata.tsv) --> summary(summary)
  sc2rf_recombinants(sc2rf_recombinants) --> summary(summary)
  usher_stats(usher_stats) --> summary(summary)
  summary(summary) --> linelist(linelist)
  issues_download(issues_download) --> linelist(linelist)
```


## Usage

```
nextflow run BCCDC-PHL/ncov2019-recombinant-screen \
  -profile conda --cache ~/.conda/envs \
  --run_dir </path/to/artic_analysis_run_dir> \
  --outdir </path/to/output_dir>
```

## Input

## Output
The output consists of five `.tsv` files:

```
output
├── <run_id>_recombinant_screen_false_positives.tsv
├── <run_id>_recombinant_screen_linelist.tsv
├── <run_id>_recombinant_screen_negatives.tsv
├── <run_id>_recombinant_screen_positives.tsv
└── <run_id>_recombinant_screen_summary.tsv
```