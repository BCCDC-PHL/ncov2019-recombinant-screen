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
include { usher_columns } from './modules/recombinant_screen.nf'
include { usher } from './modules/recombinant_screen.nf'
include { usher_stats } from './modules/recombinant_screen.nf'
include { usher_metadata } from './modules/recombinant_screen.nf'
include { usher_subtree } from './modules/recombinant_screen.nf'
include { usher_subtree_collapse } from './modules/recombinant_screen.nf'
include { summary } from './modules/recombinant_screen.nf'
include { linelist } from './modules/recombinant_screen.nf'


workflow {

  ch_run_dir            = Channel.fromPath(params.run_dir)

  ch_run_id             = Channel.of(file(params.run_dir).getName())

  if (params.ref == "NO_FILE") {
    ch_ref              = Channel.fromPath("${projectDir}/assets/reference.fasta")
  } else {
    ch_ref              = Channel.fromPath(params.ref)
  }

  if (params.breakpoints == "NO_FILE") {
    ch_breakpoints      = Channel.fromPath("${projectDir}/assets/breakpoints.tsv")
  } else {
    ch_breakpoints      = Channel.fromPath(params.breakpoints)
  }

  if (params.primer_bed == "NO_FILE") {
    ch_primer_bed       = Channel.fromPath("${projectDir}/assets/BCCDC-PHL_1200_v4.scheme.bed")
  } else {
    ch_primer_bed       = Channel.fromPath(params.primer_bed)
  }

  if (params.problematic_sites == "NO_FILE") {
    ch_problematic_sites  = Channel.fromPath("${projectDir}/assets/problematic_sites.vcf")
  } else {
    ch_problematic_sites  = Channel.fromPath(params.problematic_sites)
  }

  if (params.metadata == "NO_FILE") {
    ch_metadata           = Channel.fromPath(params.run_dir + "/" + "metadata.tsv")
  } else {
    ch_metadata           = Channel.fromPath(params.metadata)
  }

  ch_artic_analysis_dir = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir)


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

    usher_columns(ch_run_id.combine(usher_download.out.metadata))

    usher(fasta_to_vcf.out.combine(usher_download.out.pb))

    usher_stats(usher.out.clades.join(usher.out.placement_stats).join(issues_download.out.issue_to_lineage))

    usher_metadata(nextclade.out.metadata.join(sc2rf_recombinants.out.stats).join(usher_stats.out.clades).join(usher_stats.out.placement_stats).join(issues_download.out.issue_to_lineage))

    usher_subtree(usher.out.pb.join(usher_stats.out.strains).join(usher_metadata.out.decimal_date))

    usher_subtree_collapse(usher_subtree.out.subtrees_dir)

    summary(nextclade.out.metadata.join(sc2rf_recombinants.out.stats).join(usher_stats.out.clades).join(usher_stats.out.placement_stats))

    // linelist(summary.out.join(issues_download.out.issues))
}
