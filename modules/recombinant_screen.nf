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
    tuple val(run_id), path(sequences), path(dataset), path(metadata)

  output:
    tuple val(run_id), path("${run_id}.aln.fa"), emit: alignment
    tuple val(run_id), path("${run_id}_nextclade_qc.tsv"), emit: qc
    tuple val(run_id), path("${run_id}_nextclade_metadata.tsv"), emit: metadata

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

  # Merge QC output with metadata
  csvtk rename -t -f "seqName" -n "strain" ${run_id}_nextclade_qc.tsv \
      | csvtk merge -t -f "strain" <(csvtk rename -t -f "sample" -n "strain" ${metadata}) - \
      > ${run_id}_nextclade_metadata.tsv
  """
}

process nextclade_recombinants {

  tag { run_id }

  executor 'local'

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

  input:
    tuple val(run_id), path(sc2rf_ansi), path(sc2rf_csv), path(alignment), path(nextclade_qc_tsv), path(issues)

  output:
    tuple val(run_id), path("sc2rf_postprocess_output/${run_id}_sc2rf.fasta"), emit: fasta
    tuple val(run_id), path("sc2rf_postprocess_output/${run_id}_sc2rf.tsv"), emit: stats

  script:
  motifs = params.sc2rf_motifs == "" ? "" : "--motifs ${params.sc2rf_motifs}"
  """
  mkdir sc2rf_postprocess_output
  sc2rf_postprocess.py \
    --csv ${sc2rf_csv} \
    --ansi ${sc2rf_ansi} \
    --prefix ${run_id}_sc2rf \
    --min-len ${params.sc2rf_min_len} \
    --max-parents ${params.sc2rf_max_parents} \
    --max-breakpoints ${params.sc2rf_max_breakpoints} \
    --outdir sc2rf_postprocess_output \
    --aligned ${alignment} \
    --issues ${issues} \
    ${motifs}
  """
}

process nextclade_exclude_false_positives {

  tag { run_id }

  input:
    tuple val(run_id), path(alignment), path(sc2rf_stats)

  output:
    tuple val(run_id), path("${run_id}_recombinants_excluding_false_positives.aln.fa")

  script:
  """
  csvtk grep -t -f "sc2rf_status" -p "false_positive" -v ${sc2rf_stats} \
      | csvtk cut -t -f "strain" \
      | tail -n+2  \
      | seqkit grep -f - ${alignment} > ${run_id}_recombinants_excluding_false_positives.aln.fa
  """
}


process fasta_to_vcf {

  tag { run_id }

  input:
    tuple val(run_id), path(sc2rf_recombinants_fasta), path(problematic_sites), path(reference)

  output:
    tuple val(run_id), path("${run_id}_recombinants.vcf.gz")

  script:
  """
  cat ${reference} ${sc2rf_recombinants_fasta} > ${run_id}_sc2rf_recombinants_with_reference.fasta
  faToVcf \
    -ambiguousToN \
    -maskSites=${problematic_sites} \
    -ref='${params.nextclade_custom_ref}' \
    ${run_id}_sc2rf_recombinants_with_reference.fasta \
    ${run_id}_recombinants.vcf
  gzip -f ${run_id}_recombinants.vcf
  """
}

process usher_download {

  executor 'local'

  input:
    

  output:
    path("usher.pb.gz"), emit: pb
    path("metadata.tsv.gz"), emit: metadata
    path("version.txt"), emit: ver

  script:
  """
  wget -q -O usher.pb.gz ${params.usher_pb_url}
  wget -q -O metadata.tsv.gz ${params.usher_metadata_url}
  wget -q -O version.txt ${params.usher_version_url}
  """
}


process usher_columns {

  tag { run_id }

  input:
    tuple val(run_id), path(usher_metadata)

  output:
    tuple val(run_id), path("usher_metadata.tsv"), emit: usher_metadata
    tuple val(run_id), path("gisaid_id_map.tsv"), emit: gisaid_id_map

  script:
  gisaid_regex = "^.*EPI_ISL_[0-9|-]*"
  """
  set +o pipefail

  # Extract strains with embedded GISID IDs
  csvtk cut -t -f "strain" ${usher_metadata} \
    | grep -o -E "${gisaid_regex}" \
    | awk -F "|" '{{print \$0"\\\t"\$2}}' \
    > gisaid_id_map.tsv

  # Extract strains without GISAID IDs
  csvtk cut -t -f "strain" ${usher_metadata} \
    | tail -n+2 \
    | grep -v -E "${gisaid_regex}" \
    | awk '{{print \$0"\\\tNone"}}' \
    >> gisaid_id_map.tsv

  csvtk mutate2 -t -n "dataset" -e '"{wildcards.input}"' ${usher_metadata} \
    | csvtk rename -t -f "pangolin_lineage" -n "pango_lineage" \
    | csvtk mutate -t -f "strain" -n "gisaid_epi_isl" \
    | csvtk replace -t -f "gisaid_epi_isl" -p "(.*)" -k gisaid_id_map.tsv -r "{{kv}}" \
    > usher_metadata.tsv
  """
}

process usher {

  tag { run_id }

  input:
    tuple val(run_id), path(recombinants_vcf), path(usher_tree)

  output:
    tuple val(run_id), path("${run_id}_usher_clades.tsv"), emit: clades
    tuple val(run_id), path("${run_id}_usher_placement_stats.tsv"), emit: placement_stats
    tuple val(run_id), path("${run_id}_usher.pb.gz"), emit: pb

  script:
  """
  usher \
    -i ${usher_tree} \
    -o ${run_id}_usher.pb.gz \
    -v ${recombinants_vcf} \
    --threads ${task.cpus} \
    --outdir usher_output
  cp usher_output/clades.txt ${run_id}_usher_clades.tsv
  cp usher_output/placement_stats.tsv ${run_id}_usher_placement_stats.tsv
  """
}

process usher_stats {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(usher_clades), path(usher_placement_stats), path(issue_to_lineage)

  output:
    tuple val(run_id), path("${run_id}_usher_clades_with_header.tsv"), emit: clades
    tuple val(run_id), path("${run_id}_placement_stats_with_header.tsv"), emit: placement_stats
    tuple val(run_id), path("${run_id}_usher_strains.tsv"), emit: strains


  script:
  """
  echo -e "strain\\tusher_best_set_difference\\tusher_num_best" > ${run_id}_placement_stats_with_header.tsv
  sed 's/[[:space:]]*\$//' ${usher_placement_stats} >> ${run_id}_placement_stats_with_header.tsv

  echo -e "strain\\tusher_clade\\tusher_pango_lineage" > ${run_id}_usher_clades_with_header.tsv
  cat ${usher_clades} >> ${run_id}_usher_clades_with_header.tsv

  csvtk mutate -t -f "usher_pango_lineage" -n "usher_pango_lineage_map" ${run_id}_usher_clades_with_header.tsv \
    | csvtk replace -t -f "usher_pango_lineage_map" -p "proposed([0-9]+)" -k ${issue_to_lineage} -r "{{kv}}" \
    > ${run_id}_usher_clades_with_header.tsv.tmp
  mv ${run_id}_usher_clades_with_header.tsv.tmp ${run_id}_usher_clades_with_header.tsv

  csvtk cut -t -f "strain" ${run_id}_usher_clades_with_header.tsv | tail -n+2 > ${run_id}_usher_strains.tsv
  """
}

process usher_metadata {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(nextclade_metadata), path(sc2rf_recombinants_stats), path(usher_clades), path(usher_placement_stats), path(issue_to_lineage)

  output:
    tuple val(run_id), path("${run_id}_usher_metadata.tsv"), emit: metadata
    tuple val(run_id), path("${run_id}_usher_metadata_decimal_date.tsv"), emit: decimal_date

  script:
  default_cols = "strain,date"
  extract_cols = "clade,usher_clade,Nextclade_pango,usher_pango_lineage,sc2rf_parents,sc2rf_regions,sc2rf_breakpoints,usher_num_best"
  rename_cols = "Nextstrain_clade,Nextstrain_clade_usher,pangolin_lineage,pango_lineage_usher,parents,parents_regions,breakpoints,usher_placements"
  rename_cols_final = "clade_nextclade,clade_usher,lineage_nextclade,lineage_usher,parents,parents_regions,breakpoints,usher_placements"
  """
  csvtk merge -t -f "strain" ${nextclade_metadata} ${sc2rf_recombinants_stats} ${usher_clades} ${usher_placement_stats} \
    | csvtk cut -t -f "${default_cols},${extract_cols}" \
    | csvtk rename -t -f "${extract_cols}" -n "${rename_cols}" \
    | csvtk cut -t -f "${default_cols},${rename_cols}" \
    | csvtk rename -t -f "${rename_cols}" -n "${rename_cols_final}" \
    | csvtk replace -t -f "lineage_usher" -p "proposed([0-9]+)" -k ${issue_to_lineage} -r "{kv}" \
    > ${run_id}_usher_metadata.tsv

  date_to_decimal.py ${run_id}_usher_metadata.tsv ${run_id}_usher_metadata_decimal_date.tsv
  """
}


process usher_subtree {

  tag { run_id }

  input:
    tuple val(run_id), path(usher_tree), path(usher_stats_strains), path(usher_metadata_decimal_date)

  output:
    tuple val(run_id), path("${run_id}_subtrees"), emit: subtrees_dir

  script:
  """
  mkdir ${run_id}_subtrees
  cd ${run_id}_subtrees
  matUtils extract \
    -i ../${usher_tree} \
    --nearest-k-batch ../${usher_stats_strains}:${params.usher_subtree_k} \
    -M ../${usher_metadata_decimal_date} \
    --threads ${task.cpus}
  """
}

process usher_subtree_collapse {

  tag { run_id }

  input:
    tuple val(run_id), path(subtrees_dir)

  output:
    tuple val(run_id), path("${run_id}_subtrees_collapsed"), emit: collapse_dir

  script:
  """
  usher_collapse.py --indir ${subtrees_dir} --outdir ${run_id}_subtrees_collapsed
  """
}

process summary {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "${run_id}_recombinant_screen_summary.tsv", mode: 'copy'

  input:
    tuple val(run_id), path(nextclade_metadata), path(sc2rf_stats), path(usher_clades), path(usher_placements)

  output:
    tuple val(run_id), path("${run_id}_recombinant_screen_summary.tsv")

  script:
  """
  csvtk cut -t -f "strain,date,clade,Nextclade_pango" ${nextclade_metadata} \
    | csvtk rename -t -f "clade" -n "Nextclade_clade" \
    | csvtk merge -t --na "NA" -f "strain" - ${sc2rf_stats} \
    | csvtk merge -t -k --na "NA" -f "strain" ${usher_placements} - \
    | csvtk merge -t -k --na "NA" -f "strain" ${usher_clades} - \
    | csvtk sort -t -k "Nextclade_clade:r" \
    > ${run_id}_recombinant_screen_summary.tsv
  """
}

process linelist {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "${run_id}_recombinant_screen_*.tsv", mode: 'copy'

  input:
    tuple val(run_id), path(recombinant_screen_summary), path(issues)	

  output:
    tuple val(run_id), path("${run_id}_recombinant_screen_linelist.tsv"), emit: linelist 
    tuple val(run_id), path("${run_id}_recombinant_screen_positives.tsv"), emit: positives
    tuple val(run_id), path("${run_id}_recombinant_screen_false_positives.tsv"), emit: false_positives
    tuple val(run_id), path("${run_id}_recombinant_screen_negatives.tsv"), emit: negatives

  script:
  """
  linelist.py --input ${recombinant_screen_summary} --issues ${issues} --max-placements ${params.usher_max_placements} --outdir . 
  mv linelist.tsv ${run_id}_recombinant_screen_linelist.tsv
  mv positives.tsv ${run_id}_recombinant_screen_positives.tsv
  mv false_positives.tsv ${run_id}_recombinant_screen_false_positives.tsv
  mv negatives.tsv ${run_id}_recombinant_screen_negatives.tsv
  """
}

process report {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "", mode: 'copy'

  input:
    tuple val(run_id), path(nextclade_metadata), path(sc2rf_stats), path(usher_clades), path(usher_placements)

  output:
    tuple val(run_id), path("${run_id}_recombinant_screen_summary.tsv")

  script:
  """
  """
}