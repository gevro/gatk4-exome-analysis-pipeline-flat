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

import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/Alignment.wdl" as Alignment
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/SplitLargeReadGroup.wdl" as SplitRG
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/BamProcessing.wdl" as Processing
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {
  input {
    String base_file_name
    String final_gvcf_base_name
    Array[File] flowcell_unmapped_bams
    String sample_name
    String unmapped_bam_suffix

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
    File ref_alt
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

    String cross_check_fingerprints_by
    File? haplotype_database_file
    Float lod_threshold
    String recalibrated_bam_basename
  }

  #Float cutoff_for_large_rg_in_gb = 20.0
  Float cutoff_for_large_rg_in_gb = 0.1
  
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  Int compression_level = 2

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call Alignment.GetBwaVersion

  # Get the size of the standard reference files as well as the additional reference files needed for BWA

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in flowcell_unmapped_bams) {

    Float unmapped_bam_size = size(unmapped_bam, "GiB")

    String unmapped_bam_basename = basename(unmapped_bam, unmapped_bam_suffix)

    # QC the unmapped BAM
    call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = unmapped_bam_basename + ".unmapped.quality_yield_metrics",
        preemptible_tries = papi_settings.preemptible_tries
    }

    if (unmapped_bam_size > cutoff_for_large_rg_in_gb) {
      # Split bam into multiple smaller bams,
      # map reads to reference and recombine into one bam
      call SplitRG.SplitLargeReadGroup as SplitRG {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          bwa_version = GetBwaVersion.bwa_version,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_alt = ref_alt,
          ref_sa = ref_sa,
          ref_amb = ref_amb,
          ref_bwt = ref_bwt,
          ref_ann = ref_ann,
          ref_pac = ref_pac,
          compression_level = compression_level,
          preemptible_tries = papi_settings.preemptible_tries
      }
    }

    if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) {
      # Map reads to reference
      call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_alt = ref_alt,
          ref_sa = ref_sa,
          ref_amb = ref_amb,
          ref_bwt = ref_bwt,
          ref_ann = ref_ann,
          ref_pac = ref_pac,
          bwa_version = GetBwaVersion.bwa_version,
          compression_level = compression_level,
          preemptible_tries = papi_settings.preemptible_tries
      }
    }

    File output_aligned_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SplitRG.aligned_bam])
    
    String markilluminaadapters_metrics = select_first([SamToFastqAndBwaMemAndMba.markilluminaadapters_metrics, SplitRG.markilluminaadapters_metrics])
    
    Float mapped_bam_size = size(output_aligned_bam, "GiB")

    # QC the aligned but unsorted readgroup BAM
    # no reference as the input here is unsorted, providing a reference would cause an error
    call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_markilluminaadapters_metrics = markilluminaadapters_metrics,
        input_bam = output_aligned_bam,
        output_bam_prefix = unmapped_bam_basename + ".readgroup",
        preemptible_tries = papi_settings.preemptible_tries
    }
  }

  # Sum the read group bam sizes to approximate the aggregated bam size
  call Utils.SumFloats as SumFloats {
    input:
      sizes = mapped_bam_size,
      preemptible_tries = papi_settings.preemptible_tries
  }

  # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
  Float gb_size_cutoff_for_preemptibles = 110.0
  Boolean data_too_large_for_preemptibles = SumFloats.total_size > gb_size_cutoff_for_preemptibles

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = output_aligned_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      total_input_size = SumFloats.total_size,
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  # Sort aggregated+deduped BAM file and fix tags
  call Processing.SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  Float agg_bam_size = size(SortSampleBam.output_bam, "GiB")

  if (defined(haplotype_database_file)) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ SortSampleBam.output_bam ],
        input_bam_indexes = [SortSampleBam.output_bam_index],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = base_file_name + ".crosscheck",
        total_input_size = agg_bam_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Create list of sequences for scatter-gather parallelization
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      preemptible_tries = papi_settings.preemptible_tries
  }

  # Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      input_bam = SortSampleBam.output_bam,
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = base_file_name + ".preBqsr",
      preemptible_tries = papi_settings.agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call Processing.BaseRecalibrator as BaseRecalibrator {
      input:
        input_bam = SortSampleBam.output_bam,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        bqsr_scatter = bqsr_divisor,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call Processing.GatherBqsrReports as GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      preemptible_tries = papi_settings.preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call Processing.ApplyBQSR as ApplyBQSR {
      input:
        input_bam = SortSampleBam.output_bam,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        bqsr_scatter = bqsr_divisor,
        compression_level = compression_level,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call Processing.GatherSortedBamFiles as GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      total_input_size = agg_bam_size,
      compression_level = compression_level,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics
    Array[String] unsorted_read_group_markilluminaadapters_metrics = CollectUnsortedReadgroupBamQualityMetrics.markilluminaadapters_metrics
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_bam = GatherBamFiles.output_bam
    File output_bam_index = GatherBamFiles.output_bam_index
  }
}
