version 1.0
## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data processing according to the GATK Best Practices (June 2016)
## for human whole-genome and exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow AggregatedBamQC {
input {
    File base_recalibrated_bam
    File base_recalibrated_bam_index
    String base_name
    String sample_name
    String recalibrated_bam_base_name
    File? haplotype_database_file
    
    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File calling_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File? ref_alt
    File ref_sa
    File ref_amb
    File ref_bwt
    File ref_ann
    File ref_pac
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File dbsnp_vcf
    File dbsnp_vcf_index
    File evaluation_interval_list
    
    PapiSettings papi_settings
  }

  # QC the final BAM (consolidated after scattered BQSR)
  call QC.CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics {
    input:
      input_bam = base_recalibrated_bam,
      input_bam_index = base_recalibrated_bam_index,
      output_bam_prefix = base_name + ".readgroup",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the final BAM some more (no such thing as too much QC)
  call QC.CollectAggregationMetrics as CollectAggregationMetrics {
    input:
      input_bam = base_recalibrated_bam,
      input_bam_index = base_recalibrated_bam_index,
      output_bam_prefix = base_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  if (defined(haplotype_database_file) && defined(fingerprint_genotypes_file)) {
    # Check the sample BAM fingerprint against the sample array
    call QC.CheckFingerprint as CheckFingerprint {
      input:
        input_bam = base_recalibrated_bam,
        input_bam_index = base_recalibrated_bam_index,
        haplotype_database_file = haplotype_database_file,
        genotypes = fingerprint_genotypes_file,
        genotypes_index = fingerprint_genotypes_index,
        output_basename = base_name,
        sample = sample_name,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Generate a checksum per readgroup in the final BAM
  call QC.CalculateReadGroupChecksum as CalculateReadGroupChecksum {
    input:
      input_bam = base_recalibrated_bam,
      input_bam_index = base_recalibrated_bam_index,
      read_group_md5_filename = recalibrated_bam_base_name + ".bam.read_group_md5",
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  output {
    File read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
    File read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
    File read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics

    File calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file

    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics
    File agg_error_summary_metrics = CollectAggregationMetrics.error_summary_metrics

    File? fingerprint_summary_metrics = CheckFingerprint.summary_metrics
    File? fingerprint_detail_metrics = CheckFingerprint.detail_metrics
  }
}
