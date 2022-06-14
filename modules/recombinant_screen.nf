process issues_download {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(breakpoints)

  output:
    tuple val(run_id), path("issues.tsv"), path("issue_to_lineage.tsv")

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
  """
}

process nextclade_recombinants {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "${run_id}_nextclade_recombinants.tsv", mode: 'copy'

  input:
    tuple val(run_id), path(alignment), path(nextclade_qc)

  output:
    tuple val(run_id), path("${run_id}_nextclade_recombinants.tsv")

  script:
  """
  # Extract recombinants
  csvtk grep -t -f "clade,Nextclade_pango,qc.overallStatus" -r -p ".*(recombinant|X|bad|mediocre).*" ${nextclade_qc} \
    | csvtk grep -t -v -f "clade" -r -p ".*({params.exclude_clades}).*" \
    | csvtk cut -t -f "seqName" \
    | tail -n+2 \
    > ${run_id}_nextclade_recombinants.tsv
  # Filter Alignment
  seqkit grep -p "{params.custom_ref}" ${alignment} 1> aligned.fa
  seqkit grep -f "{output.strains}" ${alignment} 1>> aligned.fa
  """
}