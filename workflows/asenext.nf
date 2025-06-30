#!/usr/bin/env nextflow

/*
========================================================================================
    ASENext Workflow
========================================================================================
    Main workflow for ASENext pipeline
----------------------------------------------------------------------------------------
*/

include { INPUT_CHECK } from '../modules/local/input_check'
include { PREPARE_VCF } from '../modules/local/prepare_vcf'
include { CHROMOSOME_CHECK } from '../modules/local/chromosome_check'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { STAR_ALIGN_WASP } from '../modules/local/star_align_wasp/main'
include { FILTER_WASP_READS } from '../modules/local/filter_wasp_reads/main'
include { UMITOOLS_DEDUP } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../modules/nf-core/samtools/index/main'
include { BCFTOOLS_VIEW } from '../modules/local/bcftools_view/main'
include { BEAGLE5_BEAGLE } from '../modules/local/beagle5_beagle/main'
include { TABIX_INDEX } from '../modules/local/tabix_index/main'
include { PHASER } from '../modules/local/phaser/main'
include { PHASER_GENE_AE } from '../modules/local/phaser_gene_ae/main'
include { EXTRACT_ASE_GENES } from '../modules/local/extract_ase_genes/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'

workflow ASENEXT {
    take:
    ch_input                // channel: path to input CSV
    ch_multiqc_config       // channel: /path/to/multiqc_config.yaml
    ch_multiqc_custom_config // channel: /path/to/multiqc_custom_config.yaml

    main:
    ch_versions = Channel.empty()

    // Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Create explicitly structured channels for different processes
    INPUT_CHECK.out.csv.splitCsv(header:true, sep:',')
        .map { row ->
            def meta = [id: row.sample, single_end: false]
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = file(row.fastq_2)
            def reads = [fastq_1, fastq_2]  // Create a list of read files
            def vcf = file(row.vcf)
            return [meta, reads, vcf]  // This matches exactly what PREPARE_VCF expects
        }
        .set { ch_vcf_input }

    // Extract just the reads for FastQC
    ch_vcf_input
        .map { meta, reads, vcf -> [meta, reads] }
        .set { ch_reads }

    // Run FastQC
    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Prepare VCF files for STAR and Beagle - using the exact structure expected by the module
    PREPARE_VCF (
        ch_vcf_input,
        params.chromosome
    )
    ch_versions = ch_versions.mix(PREPARE_VCF.out.versions)

    // Create dummy channels for STAR_ALIGN_WASP
    ch_star_ignore_sjdbgtf = Channel.value(false)
    ch_seq_platform = Channel.value("illumina")
    ch_seq_center = Channel.value(false)

    // Create properly formatted index and gtf channels
    ch_star_index = Channel.fromPath(params.star_index).map { file -> [[:], file] }
    ch_gtf = Channel.fromPath(params.gtf).map { file -> [[:], file] }

    // Combine reads with VCF for STAR alignment
    ch_reads_with_vcf = ch_reads.join(PREPARE_VCF.out.vcf_star)

    // Run STAR alignment with WASP and VCF file using our custom module
    STAR_ALIGN_WASP (
        ch_reads_with_vcf.map { meta, reads, vcf ->
            // Return just meta and reads for the first input
            [meta, reads]
        },
        ch_star_index,             // index (with meta)
        ch_gtf,                    // gtf (with meta)
        ch_reads_with_vcf.map { meta, reads, vcf -> vcf }, // Extract just the VCF file
        ch_star_ignore_sjdbgtf,    // star_ignore_sjdbgtf
        ch_seq_platform,           // seq_platform
        ch_seq_center              // seq_center
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_WASP.out.versions)

    // First, filter WASP reads BEFORE UMI deduplication
    FILTER_WASP_READS (
        STAR_ALIGN_WASP.out.bam_sorted_aligned.map { meta, bam -> [ meta, bam, [] ] }  // Add empty bai
    )
    ch_versions = ch_versions.mix(FILTER_WASP_READS.out.versions)

    // Sort the filtered BAM files first
    ch_dummy_fasta = Channel.value([[:], []])  // Empty fasta reference

    SAMTOOLS_SORT (
        FILTER_WASP_READS.out.bam,
        ch_dummy_fasta  // empty fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // Index the sorted BAM
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Now run UMI deduplication on sorted, indexed BAM
    // Combine sorted BAM with its index for UMITOOLS_DEDUP
    ch_bam_with_index = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)

    UMITOOLS_DEDUP (
        ch_bam_with_index,       // bam with bai
        Channel.value(false)     // get_output_stats
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

    // Index the final deduplicated BAM
    SAMTOOLS_INDEX_FINAL (
        UMITOOLS_DEDUP.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FINAL.out.versions)

    // Extract chromosome from VCF using our custom BCFTOOLS_VIEW module
    BCFTOOLS_VIEW (
        PREPARE_VCF.out.vcf_beagle,
        params.chromosome
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    // Format beagle reference and map files
    ch_beagle_ref = params.beagle_ref ? Channel.fromPath(params.beagle_ref) : Channel.empty()
    ch_beagle_map = params.beagle_map ? Channel.fromPath(params.beagle_map) : Channel.empty()

    // Phase with Beagle
    BEAGLE5_BEAGLE (
        BCFTOOLS_VIEW.out.vcf,                  // vcf
        ch_beagle_ref,                          // refpanel
        ch_beagle_map                           // genmap
    )
    ch_versions = ch_versions.mix(BEAGLE5_BEAGLE.out.versions)

    // Index the Beagle output VCF
    TABIX_INDEX (
        BEAGLE5_BEAGLE.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_INDEX.out.versions)

    // Combine final BAM with its index for phaser
    ch_final_bam_with_index = UMITOOLS_DEDUP.out.bam.join(SAMTOOLS_INDEX_FINAL.out.bai)

    // Run phaser
    PHASER (
        TABIX_INDEX.out.vcf_indexed,
        ch_final_bam_with_index
    )
    ch_versions = ch_versions.mix(PHASER.out.versions)

    // Format gene_features
    ch_gene_features = Channel.fromPath(params.gene_features)

    // Run phaser_gene_ae
    PHASER_GENE_AE (
        PHASER.out.counts,
        ch_gene_features
    )
    ch_versions = ch_versions.mix(PHASER_GENE_AE.out.versions)

    // Extract ASE genes
    EXTRACT_ASE_GENES (
        PHASER_GENE_AE.out.ae
    )
    ch_versions = ch_versions.mix(EXTRACT_ASE_GENES.out.versions)

    // Collect all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // Run MultiQC - simplified to avoid dummy file issues
    MULTIQC (
        ch_multiqc_files.collect(),           // multiqc_files
        ch_multiqc_config,                    // multiqc_config
        ch_multiqc_custom_config,             // extra_multiqc_config
        [],                                   // multiqc_logo (empty)
        [],                                   // replace_names (empty)
        []                                    // sample_names (empty)
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions.ifEmpty(null)
}
