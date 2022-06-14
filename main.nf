#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { issues_download } from './modules/recombinant_screen.nf'
include { nextclade_dataset } from './modules/recombinant_screen.nf'
include { identify_complete_genomes } from './modules/recombinant_screen.nf'
include { prepare_multi_fasta } from './modules/recombinant_screen.nf'
include { nextclade } from './modules/recombinant_screen.nf'
include { nextclade_recombinants } from './modules/recombinant_screen.nf'


workflow {

  ch_run_id = Channel.of(params.run_id)
  ch_breakpoints = Channel.fromPath(params.breakpoints)
  ch_artic_analysis_dir = Channel.fromPath(params.artic_analysis_dir)

  main:
    issues_download(ch_run_id.combine(ch_breakpoints))

    nextclade_dataset(ch_run_id)

    identify_complete_genomes(ch_run_id.combine(ch_artic_analysis_dir))

    prepare_multi_fasta(identify_complete_genomes.out.combine(ch_artic_analysis_dir))

    nextclade(prepare_multi_fasta.out.join(nextclade_dataset.out))

    nextclade_recombinants(nextclade.out)
}
