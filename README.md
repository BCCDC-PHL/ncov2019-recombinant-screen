# Recombinant Screen for SARS-CoV-2

This pipeline is based directly on [ktmeaton/ncov-recombinant](https://github.com/ktmeaton/ncov-recombinant).
We've adapted the original Snakemake pipeline into a Nextflow pipeine, and made minor adjustments to fit our
systems and SARS-CoV-2 sequence datasets.

## Overview

```mermaid
flowchart TD
    p0((Channel.of))
    p1((Channel.fromPath))
    p2((Channel.fromPath))
    p3((Channel.fromPath))
    p4((Channel.fromPath))
    p5((Channel.fromPath))
    p6((Channel.fromPath))
    p7([combine])
    p8[issues_download]
    p9[nextclade_dataset]
    p10([combine])
    p11[identify_complete_genomes]
    p12([combine])
    p13[prepare_multi_fasta]
    p14([join])
    p15([combine])
    p16[nextclade]
    p17([join])
    p18[nextclade_recombinants]
    p19(( ))
    p20([combine])
    p21[sc2rf]
    p22([join])
    p23([join])
    p24([join])
    p25[sc2rf_recombinants]
    p26([join])
    p27[nextclade_exclude_false_positives]
    p28(( ))
    p29([combine])
    p30([combine])
    p31[fasta_to_vcf]
    p32[usher_download]
    p33([combine])
    p34[usher]
    p35([join])
    p36[usher_stats]
    p37(( ))
    p38([join])
    p39([join])
    p40[summary]
    p41([join])
    p42[linelist]
    p43(( ))
    p44(( ))
    p45(( ))
    p46(( ))
    p0 -->|ch_run_id| p7
    p1 -->|ch_ref| p30
    p2 -->|ch_breakpoints| p7
    p3 -->|ch_primer_bed| p20
    p4 -->|ch_problematic_sites| p29
    p5 -->|ch_artic_analysis_dir| p10
    p6 -->|ch_metadata| p15
    p7 --> p8
    p8 --> p24
    p8 --> p35
    p0 -->|ch_run_id| p9
    p9 --> p14
    p0 -->|ch_run_id| p10
    p10 --> p11
    p11 --> p12
    p5 -->|ch_artic_analysis_dir| p12
    p12 --> p13
    p13 --> p14
    p14 --> p15
    p15 --> p16
    p16 --> p17
    p16 --> p17
    p16 --> p38
    p17 --> p18
    p18 --> p20
    p18 --> p19
    p20 --> p21
    p21 --> p23
    p16 --> p22
    p16 --> p22
    p22 --> p23
    p23 --> p24
    p24 --> p25
    p25 --> p29
    p25 --> p26
    p16 --> p26
    p26 --> p27
    p27 --> p28
    p29 --> p30
    p30 --> p31
    p31 --> p33
    p32 --> p33
    p33 --> p34
    p34 --> p35
    p35 --> p36
    p36 --> p39
    p36 --> p37
    p25 --> p38
    p38 --> p39
    p39 --> p40
    p40 --> p41
    p8 --> p41
    p41 --> p42
    p42 --> p46
    p42 --> p45
    p42 --> p44
    p42 --> p43
```


## Usage

```
nextflow run BCCDC-PHL/ncov2019-recombinant-screen \
  
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