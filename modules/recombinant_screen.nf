process issues_download {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(breakpoints)

  output:
    tuple val(run_id), path("issues.tsv"), emit: issues
    tuple val(run_id), path("issue_to_lineage.tsv"), emit: issue_to_lineage

  script:
  """
  issues.py --breakpoints ${breakpoints} > issues.tsv
  csvtk cut -t -f "issue,lineage" issues.tsv | tail -n+2 > issue_to_lineage.tsv
  """
}

process nextclade_dataset {

  tag { run_id }

  executor 'local'

  input:
    val(run_id)

  output:
    tuple val(run_id), path("nextclade_${params.nextclade_dataset}")

  script:
  """
  nextclade dataset get --name ${params.nextclade_dataset} --tag ${params.nextclade_tag} --output-dir nextclade_${params.nextclade_dataset}
  """
}

process identify_complete_genomes {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(artic_analysis_dir)

  output:
    tuple val(run_id), path("${run_id}_qc_pass_sample_ids.csv")

  script:
  """
  tail -qn+2 ${artic_analysis_dir}/${run_id}.qc.csv | awk -F ',' '\$3 > ${params.minimum_genome_completeness}' | cut -d ',' -f 1 | grep -v '^POS' | grep -v '^NEG' > ${run_id}_qc_pass_sample_ids.csv
  """
}

process prepare_multi_fasta {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(qc_pass_sample_ids), path(artic_analysis_dir)

  output:
    tuple val(run_id), path("${run_id}_qc_pass_seqs.fa")

  script:
  """
  while read -r sample_id; do
    cat ${artic_analysis_dir}/${params.artic_consensus_subdir}/\${sample_id}.consensus.fa >> ${run_id}_qc_pass_seqs.fa
  done < ${qc_pass_sample_ids}
  """
}

process nextclade {

  tag { run_id }

  input:
    tuple val(run_id), path(sequences), path(dataset)

  output:
    tuple val(run_id), path("${run_id}.aln.fa"), path("${run_id}_nextclade_qc.tsv")

  script:
  """
  nextclade run \
    --jobs ${task.cpus} \
    --input-fasta ${sequences} \
    --input-dataset ${dataset} \
    --include-reference \
    --output-dir ${run_id}_nextclade \
    --output-tsv ${run_id}_nextclade_qc.tsv \
    --output-fasta ${run_id}.aln.fa \
    --output-basename ${run_id} \
    > nextclade.log 2>&1

  # Rename ref strain
  seqkit replace -p "${params.nextclade_ref}" -r "${params.nextclade_custom_ref}" ${run_id}.aln.fa > ${run_id}.aln.fa.tmp
  mv ${run_id}.aln.fa.tmp ${run_id}.aln.fa
  """
}

process nextclade_recombinants {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "${run_id}_nextclade_recombinants.tsv", mode: 'copy'

  input:
    tuple val(run_id), path(alignment), path(nextclade_qc)

  output:
    tuple val(run_id), path("${run_id}_recombinants.aln.fa"), emit: alignment
    tuple val(run_id), path("${run_id}_nextclade_recombinants.tsv"), emit: recombinants_list

  script:
  """
  # Extract recombinants
  csvtk grep -t -f "clade,Nextclade_pango,qc.overallStatus" -r -p ".*(recombinant|X|bad|mediocre).*" ${nextclade_qc} \
    | csvtk grep -t -v -f "clade" -r -p ".*(${params.nextclade_exclude_clades}).*" \
    | csvtk cut -t -f "seqName" \
    | tail -n+2 \
    > ${run_id}_nextclade_recombinants.tsv
  # Filter Alignment
  seqkit grep -p "${params.nextclade_custom_ref}" ${alignment} 1> aligned.fa
  seqkit grep -f ${run_id}_nextclade_recombinants.tsv ${alignment} 1>> ${run_id}_recombinants.aln.fa
  """
}

process sc2rf {

  tag { run_id }

  input:
    tuple val(run_id), path(alignment), path(primer_bed)

  output:
    tuple val(run_id), path("${run_id}_sc2rf_ansi.txt"), path("${run_id}_sc2rf.csv")

  script:
  sc2rf_args = "--ansi --parents 2-4 --breakpoints 0-4 --unique 1 --max-ambiguous 20 --max-intermission-length 3 --max-intermission-count 3"
  """
  sc2rf \
    ${sc2rf_args} \
    --primers ${primer_bed} \
    --max-name-length ${params.sc2rf_max_name_length} \
    --clades ${params.sc2rf_clades} \
    --mutation-threshold ${params.sc2rf_mutation_threshold} \
    --csvfile ${run_id}_sc2rf.csv \
    ${alignment} \
    > ${run_id}_sc2rf_ansi.txt
  """
}

process sc2rf_recombinants {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "sc2rf_postprocess_output/${run_id}_sc2rf.*", mode: 'copy', saveAs: { filename -> filename.split("/").last() }

  input:
    tuple val(run_id), path(sc2rf_ansi), path(sc2rf_csv), path(alignment), path(nextclade_qc_tsv), path(issues)

  output:
    tuple val(run_id), path("sc2rf_postprocess_output/${run_id}_sc2rf.fasta"), emit: fasta
    tuple val(run_id), path("sc2rf_postprocess_output/${run_id}_sc2rf.tsv"), emit: tsv

  script:
  motifs = params.sc2rf_motifs == "" ? "" : "--motifs ${params.sc2rf_motifs}"
  """
  mkdir sc2rf_postprocess_output
  sc2rf_postprocess.py \
    --prefix ${run_id}_sc2rf \
    --csv ${sc2rf_csv} \
    --ansi ${sc2rf_ansi} \
    --min-len ${params.sc2rf_min_len} \
    --max-parents ${params.sc2rf_max_parents} \
    --max-breakpoints ${params.sc2rf_max_breakpoints} \
    --outdir sc2rf_postprocess_output \
    --aligned ${alignment} \
    --custom-ref "${params.nextclade_custom_ref}" \
    --issues ${issues} \
    --qc ${nextclade_qc_tsv} \
    ${motifs}
  """
}

process fasta_to_vcf {

  tag { run_id }

  input:
    tuple val(run_id), path(sc2rf_recombinants_fasta), path(problematic_sites)

  output:
    tuple val(run_id), path("${run_id}_recombinants.vcf.gz")

  script:
  """
  faToVcf \
    -ambiguousToN \
    -maskSites=${problematic_sites} \
    -ref='${params.nextclade_custom_ref}' \
    ${sc2rf_recombinants_fasta} \
    ${run_id}_recombinants.vcf
  gzip -f ${run_id}_recombinants.vcf
  """
}

process usher_download {

  tag {}

  input:
    

  output:
    

  script:
  """
  wget -q -O usher.pb.gz ${params.usher_pb_url}
  wget -q -O metadata.tsv.gz ${params.usher_metadata_url}
  wget -q -O version.txt ${params.usher_version_url}
  # csvtk mutate2 -t -n "dataset" -e '"{wildcards.input}"' {output.metadata}.gz 1> {output.metadata} 2>> {log};
  # rm -f {output.metadata}.gz
  """
}