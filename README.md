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
  nextclade(nextclade) -- qc --> nextclade_recombinants(nextclade_recombinants)
  nextclade(nextclade) -- alignment --> nextclade_recombinants(nextclade_recombinants)
  nextclade_recombinants(nextclade_recombinants) -- alignment --> sc2rf(sc2rf)
  primers --> sc2rf(sc2rf)
  sc2rf(sc2rf) -- ansi --> sc2rf_recombinants(sc2rf_recombinants)
  sc2rf(sc2rf) -- csv --> sc2rf_recombinants(sc2rf_recombinants)
  nextclade(nextclade) -- alignment --> sc2rf(sc2rf)
  issues_download(issues_download) -- issues --> sc2rf(sc2rf)
  nextclade(nextclade) --> nextclade_exclude_false_positives(nextclade_exclude_false_positives)
  sc2rf_recombinants(sc2rf_recombinants) -- stats --> nextclade_exclude_false_positives(nextclade_exclude_false_positives)
  sc2rf_recombinants(sc2rf_recombinants) -- fasta --> fasta_to_vcf(fasta_to_vcf)
  problematic_sites(problematic_sites.vcf) --> fasta_to_vcf(fasta_to_vcf)
  ref(ref.fa) --> fasta_to_vcf(fasta_to_vcf)
  usher_download(usher_download) -- pb_tree --> usher(usher)
  usher_download(usher_download) -- metadata --> usher_columns(usher_columns)
  fasta_to_vcf(fasta_to_vcf) -- recombinants_vcf --> usher(usher)
  usher(usher) -- clades --> usher_stats(usher_stats)
  usher(usher) -- placement_stats --> usher_stats(usher_stats)
  nextclade(nextclade) -- metadata --> usher_metadata(usher_metadata)
  sc2rf_recombinants(sc2rf_recombinants) -- stats --> usher_metadata(usher_metadata)
  usher_stats(usher_stats) -- clades --> usher_metadata(usher_metadata)
  nextclade(nextclade) -- metadata --> usher_metadata(usher_metadata)
  sc2rf_recombinants(sc2rf_recombinants) -- stats --> usher_metadata(usher_metadata)
  usher_stats(usher_stats) -- clades --> usher_metadata(usher_metadata)
  usher_stats(usher_stats) -- placement_stats --> usher_metadata(usher_metadata)
  issues_download(issues_download) -- issue_to_lineage --> usher_metadata(usher_metadata)
  usher(usher) -- pb_tree --> usher_subtree(usher_subtree)
  usher_stats(usher_stats) -- strains --> usher_subtree(usher_subtree)
  usher_metadata(usher_metadata) -- decimal_to_date --> usher_subtree(usher_subtree)
  issues_download(issues_download) -- issue_to_lineage--> usher_stats(usher_stats)
  nextclade(nextclade) -- metadata --> summary(summary)
  sc2rf_recombinants(sc2rf_recombinants) --> summary(summary)
  usher_stats(usher_stats) --> summary(summary)
  summary(summary) --> linelist(linelist)
  issues_download(issues_download) --> linelist(linelist)
  geo_resolutions(geo_resolutions.json) --> usher_subtree_collapse(usher_subtree_collapse)
  usher_subtree(usher_subtree) -- subtrees_dir --> usher_subtree_collapse(usher_subtree_collapse)
  summary(summary) -- summary --> linelist(linelist)
  issues_download(issues_download) -- issues --> linelist(linelist)
  linelist(linelist) --> linelist_tsv(linelist.tsv)
  linelist(linelist) --> positives_tsv(positives.tsv)
  linelist(linelist) --> negatives_tsv(negatives.tsv)
  linelist(linelist) --> false_positives_tsv(false_positives.tsv)
  style nextclade fill:#4287f5
  style usher fill:#4287f5
  style sc2rf fill:#4287f5
  style linelist_tsv fill: #53cc43
  style positives_tsv fill: #53cc43
  style negatives_tsv fill: #53cc43
  style false_positives_tsv fill: #53cc43
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
