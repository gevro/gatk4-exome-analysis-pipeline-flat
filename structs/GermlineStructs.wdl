version 1.0

struct CrossSpeciesContaminationReferences {
  File filter_bwa_image
  File kmer_file
  File meats_bwa_image
  File meats_fasta
  File meats_fasta_dict
  File meats_taxonomy_file
  File microbe_bwa_image
  File microbe_fasta
  File microbe_fasta_dict
  File microbe_taxonomy_file
  File normalization_file
  File metrics_script_file
  Float score_min_identity
  Int reads_after_downsampling
}

struct PapiSettings {
  Int preemptible_tries
  Int agg_preemptible_tries
}
