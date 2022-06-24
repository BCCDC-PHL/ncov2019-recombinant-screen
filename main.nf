#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { issues_download } from './modules/recombinant_screen.nf'
include { nextclade_dataset } from './modules/recombinant_screen.nf'
include { identify_complete_genomes } from './modules/recombinant_screen.nf'
include { prepare_multi_fasta } from './modules/recombinant_screen.nf'
include { nextclade } from './modules/recombinant_screen.nf'
include { nextclade_recombinants } from './modules/recombinant_screen.nf'
include { sc2rf } from './modules/recombinant_screen.nf'
include { sc2rf_recombinants } from './modules/recombinant_screen.nf'
include { nextclade_exclude_false_positives } from './modules/recombinant_screen.nf'
include { fasta_to_vcf } from './modules/recombinant_screen.nf'
include { usher_download } from './modules/recombinant_screen.nf'
include { usher } from './modules/recombinant_screen.nf'
include { summary } from './modules/recombinant_screen.nf'


workflow {

  ch_run_id             = Channel.of(params.run_id)
  ch_ref                = Channel.fromPath(params.ref)
  ch_breakpoints        = Channel.fromPath(params.breakpoints)
  ch_primer_bed         = Channel.fromPath(params.primer_bed)
  ch_problematic_sites  = Channel.fromPath(params.problematic_sites)
  ch_artic_analysis_dir = Channel.fromPath(params.artic_analysis_dir)
  ch_metadata           = Channel.fromPath(params.metadata)

  main:
    issues_download(ch_run_id.combine(ch_breakpoints))

    nextclade_dataset(ch_run_id)

    identify_complete_genomes(ch_run_id.combine(ch_artic_analysis_dir))

    prepare_multi_fasta(identify_complete_genomes.out.combine(ch_artic_analysis_dir))

    nextclade(prepare_multi_fasta.out.join(nextclade_dataset.out).combine(ch_metadata))

    nextclade_recombinants(nextclade.out.alignment.join(nextclade.out.qc))

    sc2rf(nextclade_recombinants.out.alignment.combine(ch_primer_bed))

    sc2rf_recombinants(sc2rf.out.join(nextclade.out.alignment.join(nextclade.out.qc)).join(issues_download.out.issues))

    nextclade_exclude_false_positives(nextclade.out.alignment.join(sc2rf_recombinants.out.stats))

    fasta_to_vcf(sc2rf_recombinants.out.fasta.combine(ch_problematic_sites).combine(ch_ref))
    
    usher_download()

    usher(fasta_to_vcf.out.combine(usher_download.out))

    summary(nextclade.out.metadata.join(sc2rf_recombinants.out.stats))
}
