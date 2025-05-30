# ASENext

## Overview
ASENext is a Nextflow pipeline for allele-specific expression (ASE) analysis using STAR-WASP for alignment, UMI-tools for deduplication, and phaser for haplotype phasing and ASE detection.

## Features
- STAR alignment with WASP mode for allele-specific mapping
- UMI-based deduplication
- Chromosome-specific analysis with configurable chromosome selection
- Beagle phasing integration
- Phaser-based allele-specific expression analysis
- Comprehensive QC with FastQC and MultiQC reporting

## Requirements
- Nextflow (>=21.10.3)
- Singularity or Docker
- Reference genome and annotation files
- Beagle reference panel and genetic map (for phasing)

## Usage

### Basic usage
```bash
nextflow run ASENext --input samples.csv --outdir results --chromosome chr11
```

### Input format
The pipeline requires a CSV file with the following columns:
```csv
sample,fastq_1,fastq_2,vcf
SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,/path/to/sample1.vcf.gz
```

### Parameters

#### Required parameters
- `--input`: Path to input CSV file
- `--outdir`: Path to output directory

#### Reference genome parameters
- `--genome`: Name of iGenomes reference (e.g., 'GRCh38')
- `--fasta`: Path to reference genome FASTA file (if not using `--genome`)
- `--gtf`: Path to GTF annotation file (if not using `--genome`)
- `--star_index`: Path to STAR index directory (if not using `--genome`)
- `--gene_features`: Path to BED file with gene features for phaser_gene_ae

#### Chromosome and phasing parameters
- `--chromosome`: Chromosome to analyze (default: 'chr11')
- `--beagle_ref`: Path to Beagle reference panel VCF
- `--beagle_map`: Path to Beagle genetic map file

#### UMI parameters
- `--umi_separator`: UMI separator character in read IDs (default: ':')

## Pipeline steps
1. Input validation and VCF preparation
2. FastQC for raw reads
3. STAR alignment with WASP mode
4. Filtering of WASP-passing reads
5. UMI-based deduplication
6. Sorting and indexing of BAM files
7. Chromosome extraction from VCF
8. Beagle phasing
9. Phaser for haplotype-level expression
10. Phaser_gene_ae for gene-level ASE
11. Extraction of ASE genes
12. MultiQC report generation

## Output
The pipeline organizes outputs by sample name in the specified output directory:
- `fastqc/`: FastQC reports
- `star/`: STAR alignment results
- `wasp/`: WASP-filtered BAM files
- `umi/`: UMI-deduplicated BAM files
- `beagle/`: Phased VCF files
- `phaser/`: Phaser results
- `ase/`: Allele-specific expression results
- `multiqc/`: MultiQC report

## Credits
- Author: Abu Saadat
- Pipeline framework: nf-core
- Tools: STAR, UMI-tools, samtools, bcftools, Beagle, phaser

## Contributing

ASENext is under active development and we welcome contributions!
If you find a bug, have an idea to improve it, or want to help implement new features (like better sex-chromosome support), feel free to open an issue or submit a pull request.

Let’s build this together.
