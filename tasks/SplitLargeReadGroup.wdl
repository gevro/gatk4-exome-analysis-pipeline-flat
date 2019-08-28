version 1.0
## Copyright Broad Institute, 2018
##
## This WDL pipeline implements a split of large readgroups for human whole-genome and exome sequencing data.
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
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/BamProcessing.wdl" as Processing
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/tasks/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/gevro/gatk4-exome-analysis-pipeline-flat/master/structs/GermlineStructs.wdl"

workflow SplitLargeReadGroup {
  input {
    File input_bam

    String bwa_commandline
    String bwa_version
    String output_bam_basename

    # ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File ref_alt
    File ref_sa
    File ref_amb
    File ref_bwt
    File ref_ann
    File ref_pac

    Int compression_level
    Int preemptible_tries
    Int reads_per_file = 48000000

  }

  call Alignment.SamSplitter as SamSplitter {
    input :
      input_bam = input_bam,
      n_reads = reads_per_file,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    Float current_unmapped_bam_size = size(unmapped_bam, "GiB")
    String current_name = basename(unmapped_bam, ".bam")

    call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        #output_bam_basename = current_name,
        output_bam_basename = basename(input_bam,".bam") + "." + current_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_alt = ref_alt,
        ref_sa = ref_sa,
        ref_amb = ref_amb,
        ref_bwt = ref_bwt,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        bwa_version = bwa_version,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries
    }

    Float current_mapped_size = size(SamToFastqAndBwaMemAndMba.output_bam, "GiB")
  }

  call Utils.SumFloats as SumSplitAlignedSizes {
    input:
      sizes = current_mapped_size,
      preemptible_tries = preemptible_tries
  }

  call Processing.GatherUnsortedBamFiles as GatherMonolithicBamFile {
    input:
      input_bams = SamToFastqAndBwaMemAndMba.output_bam,
      total_input_size = SumSplitAlignedSizes.total_size,
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }
  output {
    File aligned_bam = GatherMonolithicBamFile.output_bam
    String markilluminaadapters_metrics = "${sep=',' SamToFastqAndBwaMemAndMba.markilluminaadapters_metrics}"
  }
}
