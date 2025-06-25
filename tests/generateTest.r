#!/usr/bin/env Rscript

# ASENext Test Data Generator
# This script creates small test datasets for the ASENext pipeline

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)

# Create output directory
dir.create("test_data", showWarnings = FALSE)
dir.create("test_data/fastq", showWarnings = FALSE)
dir.create("test_data/reference", showWarnings = FALSE)
dir.create("test_data/variants", showWarnings = FALSE)

# =============================================================================
# 1. Create a mini reference genome (chromosome 11 subset)
# =============================================================================

# Create a small reference sequence for chr11 (1MB region)
set.seed(123)
chr11_seq <- sample(c("A", "T", "G", "C"), size = 1000000, replace = TRUE, 
                   prob = c(0.25, 0.25, 0.25, 0.25))
chr11_seq <- paste(chr11_seq, collapse = "")

# Write reference FASTA
ref_fasta <- DNAStringSet(chr11_seq)
names(ref_fasta) <- "chr11"
writeXStringSet(ref_fasta, "test_data/reference/chr11_mini.fa")

cat("Created mini reference genome: test_data/reference/chr11_mini.fa\n")

# =============================================================================
# 2. Create a simple GTF annotation file
# =============================================================================

# Create some test genes
genes <- data.frame(
  seqname = "chr11",
  source = "test",
  feature = "gene",
  start = c(10000, 50000, 150000, 250000, 350000),
  end = c(30000, 80000, 180000, 280000, 380000),
  score = ".",
  strand = c("+", "-", "+", "+", "-"),
  frame = ".",
  attribute = paste0('gene_id "GENE', 1:5, '"; gene_name "TEST', 1:5, '"; gene_type "protein_coding";')
)

# Create exons for each gene
gtf_lines <- c()
for(i in 1:nrow(genes)) {
  # Add gene line
  gtf_lines <- c(gtf_lines, paste(genes[i,], collapse = "\t"))
  
  # Add transcript line
  transcript_line <- genes[i,]
  transcript_line$feature <- "transcript"
  transcript_line$attribute <- paste0(genes[i, "attribute"], ' transcript_id "GENE', i, '.1";')
  gtf_lines <- c(gtf_lines, paste(transcript_line, collapse = "\t"))
  
  # Add 2-3 exons per gene
  gene_start <- genes[i, "start"]
  gene_end <- genes[i, "end"]
  gene_length <- gene_end - gene_start
  
  if(gene_length > 5000) {
    # Multi-exon gene
    exon1_start <- gene_start
    exon1_end <- gene_start + 1000
    exon2_start <- gene_end - 1000
    exon2_end <- gene_end
    
    for(exon_num in 1:2) {
      exon_line <- genes[i,]
      exon_line$feature <- "exon"
      if(exon_num == 1) {
        exon_line$start <- exon1_start
        exon_line$end <- exon1_end
      } else {
        exon_line$start <- exon2_start
        exon_line$end <- exon2_end
      }
      exon_line$attribute <- paste0(genes[i, "attribute"], 
                                   ' transcript_id "GENE', i, '.1"; exon_number "', exon_num, '";')
      gtf_lines <- c(gtf_lines, paste(exon_line, collapse = "\t"))
    }
  } else {
    # Single exon gene
    exon_line <- genes[i,]
    exon_line$feature <- "exon"
    exon_line$attribute <- paste0(genes[i, "attribute"], 
                                 ' transcript_id "GENE', i, '.1"; exon_number "1";')
    gtf_lines <- c(gtf_lines, paste(exon_line, collapse = "\t"))
  }
}

writeLines(gtf_lines, "test_data/reference/chr11_mini.gtf")
cat("Created GTF annotation: test_data/reference/chr11_mini.gtf\n")

# =============================================================================
# 3. Create BED file for gene features (needed for phaser_gene_ae)
# =============================================================================

bed_lines <- c()
for(i in 1:nrow(genes)) {
  bed_line <- paste(genes[i, "seqname"], 
                   genes[i, "start"] - 1,  # BED is 0-based
                   genes[i, "end"],
                   paste0("GENE", i),
                   "0",  # score
                   genes[i, "strand"],
                   sep = "\t")
  bed_lines <- c(bed_lines, bed_line)
}

writeLines(bed_lines, "test_data/reference/chr11_mini_genes.bed")
cat("Created BED file: test_data/reference/chr11_mini_genes.bed\n")

# =============================================================================
# 4. Create test VCF files with heterozygous SNPs
# =============================================================================

create_vcf <- function(sample_name, n_variants = 50) {
  set.seed(123 + match(sample_name, c("sample1", "sample2")))
  
  # Generate random SNP positions
  positions <- sort(sample(20000:950000, n_variants))
  
  # Create VCF header
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr11,length=1000000>",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name, sep = "\t")
  )
  
  # Generate variant lines
  vcf_lines <- c()
  for(i in 1:length(positions)) {
    pos <- positions[i]
    ref_allele <- sample(c("A", "T", "G", "C"), 1)
    alt_allele <- sample(setdiff(c("A", "T", "G", "C"), ref_allele), 1)
    
    # Create heterozygous genotype with realistic read depths
    ref_depth <- sample(15:35, 1)
    alt_depth <- sample(15:35, 1)
    total_depth <- ref_depth + alt_depth
    
    variant_line <- paste(
      "chr11",                    # CHROM
      pos,                        # POS
      ".",                        # ID
      ref_allele,                 # REF
      alt_allele,                 # ALT
      "60",                       # QUAL
      "PASS",                     # FILTER
      paste0("AC=1;AN=2"),        # INFO
      "GT:AD:DP",                 # FORMAT
      paste0("0/1:", ref_depth, ",", alt_depth, ":", total_depth),  # Sample
      sep = "\t"
    )
    vcf_lines <- c(vcf_lines, variant_line)
  }
  
  # Write VCF file
  vcf_content <- c(vcf_header, vcf_lines)
  vcf_file <- paste0("test_data/variants/", sample_name, ".vcf")
  writeLines(vcf_content, vcf_file)
  
  # Compress and index
  system(paste("bgzip -f", vcf_file))
  system(paste("tabix -f -p vcf", paste0(vcf_file, ".gz")))
  
  cat("Created VCF file:", paste0(vcf_file, ".gz"), "\n")
  return(paste0(vcf_file, ".gz"))
}

# Create VCF files for two samples
vcf_files <- sapply(c("sample1", "sample2"), create_vcf)

# =============================================================================
# 5. Create synthetic FASTQ files with UMIs
# =============================================================================

create_fastq <- function(sample_name, n_reads = 1000) {
  set.seed(123 + match(sample_name, c("sample1", "sample2")) * 10)
  
  # Generate UMI sequences (8bp)
  generate_umi <- function() {
    paste(sample(c("A", "T", "G", "C"), 8, replace = TRUE), collapse = "")
  }
  
  # Generate reads mapping to gene regions
  reads_r1 <- c()
  reads_r2 <- c()
  
  for(i in 1:n_reads) {
    # Generate UMI
    umi <- generate_umi()
    
    # Generate read sequences (simplified - just random sequences)
    read1 <- paste(sample(c("A", "T", "G", "C"), 75, replace = TRUE), collapse = "")
    read2 <- paste(sample(c("A", "T", "G", "C"), 75, replace = TRUE), collapse = "")
    
    # Quality scores (simplified)
    qual <- paste(rep("I", 75), collapse = "")
    
    # Read headers with UMI in read name
    header1 <- paste0("@", sample_name, "_", sprintf("%06d", i), ":", umi, "/1")
    header2 <- paste0("@", sample_name, "_", sprintf("%06d", i), ":", umi, "/2")
    
    reads_r1 <- c(reads_r1, header1, read1, "+", qual)
    reads_r2 <- c(reads_r2, header2, read2, "+", qual)
  }
  
  # Write FASTQ files
  fastq1 <- paste0("test_data/fastq/", sample_name, "_R1.fastq")
  fastq2 <- paste0("test_data/fastq/", sample_name, "_R2.fastq")
  
  writeLines(reads_r1, fastq1)
  writeLines(reads_r2, fastq2)
  
  # Compress FASTQ files
  system(paste("gzip -f", fastq1))
  system(paste("gzip -f", fastq2))
  
  cat("Created FASTQ files:", paste0(fastq1, ".gz"), "and", paste0(fastq2, ".gz"), "\n")
  
  return(c(paste0(fastq1, ".gz"), paste0(fastq2, ".gz")))
}

# Create FASTQ files for two samples
fastq_files <- lapply(c("sample1", "sample2"), create_fastq)

# =============================================================================
# 6. Create sample sheet
# =============================================================================

samplesheet <- data.frame(
  sample = c("sample1", "sample2"),
  fastq_1 = c(fastq_files[[1]][1], fastq_files[[2]][1]),
  fastq_2 = c(fastq_files[[1]][2], fastq_files[[2]][2]),
  vcf = vcf_files,
  stringsAsFactors = FALSE
)

write.csv(samplesheet, "test_data/samplesheet.csv", row.names = FALSE, quote = FALSE)
cat("Created samplesheet: test_data/samplesheet.csv\n")

# =============================================================================
# 7. Create minimal Beagle reference files (for testing)
# =============================================================================

# Create a minimal reference VCF for Beagle
beagle_header <- c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr11,length=1000000>",
  paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "REF_SAMPLE", sep = "\t")
)

# Add a few reference variants
set.seed(456)
ref_positions <- sort(sample(10000:900000, 20))
beagle_lines <- c()

for(pos in ref_positions) {
  ref_allele <- sample(c("A", "T", "G", "C"), 1)
  alt_allele <- sample(setdiff(c("A", "T", "G", "C"), ref_allele), 1)
  
  line <- paste(
    "chr11", pos, ".", ref_allele, alt_allele, "60", "PASS", ".", "GT", "0|1",
    sep = "\t"
  )
  beagle_lines <- c(beagle_lines, line)
}

beagle_content <- c(beagle_header, beagle_lines)
writeLines(beagle_content, "test_data/variants/beagle_ref_chr11.vcf")
system("bgzip -f test_data/variants/beagle_ref_chr11.vcf")
system("tabix -f -p vcf test_data/variants/beagle_ref_chr11.vcf.gz")

# Create minimal genetic map
map_lines <- c("pos chr cM")
for(i in 1:length(ref_positions)) {
  map_line <- paste(ref_positions[i], "11", ref_positions[i] / 1000000, sep = " ")
  map_lines <- c(map_lines, map_line)
}
writeLines(map_lines, "test_data/variants/genetic_map_chr11.txt")

cat("Created Beagle reference files\n")

# =============================================================================
# 8. Create test configuration file
# =============================================================================

test_config <- '
params {
    // Test data paths
    input                 = "test_data/samplesheet.csv"
    outdir                = "test_results"
    
    // Reference files
    fasta                 = "test_data/reference/chr11_mini.fa"
    gtf                   = "test_data/reference/chr11_mini.gtf"
    gene_features         = "test_data/reference/chr11_mini_genes.bed"
    
    // You will need to build STAR index from the reference:
    // mkdir -p test_data/reference/star_index
    // STAR --runMode genomeGenerate --genomeDir test_data/reference/star_index --genomeFastaFiles test_data/reference/chr11_mini.fa --sjdbGTFfile test_data/reference/chr11_mini.gtf --genomeSAindexNbases 10
    star_index            = "test_data/reference/star_index"
    
    // Beagle files
    beagle_ref            = "test_data/variants/beagle_ref_chr11.vcf.gz"
    beagle_map            = "test_data/variants/genetic_map_chr11.txt"
    
    // Analysis parameters
    chromosome            = "chr11"
    umi_separator         = ":"
    
    // Resources for testing
    max_cpus              = 2
    max_memory            = "4.GB"
    max_time              = "2.h"
}

process {
    withName: STAR_ALIGN_WASP {
        cpus = 2
        memory = "4.GB"
    }
}
'

writeLines(test_config, "test_data/test.config")

# =============================================================================
# 9. Create setup script
# =============================================================================

setup_script <- '#!/bin/bash

# ASENext Test Data Setup Script

echo "Setting up test data for ASENext pipeline..."

# Check if required tools are available
command -v STAR >/dev/null 2>&1 || { echo "STAR is required but not installed. Aborting." >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools is required but not installed. Aborting." >&2; exit 1; }

# Create STAR index
echo "Creating STAR index..."
mkdir -p test_data/reference/star_index
STAR --runMode genomeGenerate \
     --genomeDir test_data/reference/star_index \
     --genomeFastaFiles test_data/reference/chr11_mini.fa \
     --sjdbGTFfile test_data/reference/chr11_mini.gtf \
     --genomeSAindexNbases 10 \
     --runThreadN 2

# Create samtools index for reference
echo "Indexing reference FASTA..."
samtools faidx test_data/reference/chr11_mini.fa

echo "Test data setup complete!"
echo ""
echo "To run the test:"
echo "nextflow run main.nf -c test_data/test.config"
echo ""
echo "Or with Docker:"
echo "nextflow run main.nf -c test_data/test.config -profile docker"
'

writeLines(setup_script, "test_data/setup_test_data.sh")
system("chmod +x test_data/setup_test_data.sh")

cat("\n=============================================================================\n")
cat("Test data creation complete!\n")
cat("=============================================================================\n")
cat("Created files:\n")
cat("- test_data/reference/chr11_mini.fa (mini reference genome)\n")
cat("- test_data/reference/chr11_mini.gtf (gene annotations)\n") 
cat("- test_data/reference/chr11_mini_genes.bed (gene features for phaser)\n")
cat("- test_data/fastq/sample[1,2]_R[1,2].fastq.gz (test FASTQ files)\n")
cat("- test_data/variants/sample[1,2].vcf.gz (test VCF files)\n")
cat("- test_data/variants/beagle_ref_chr11.vcf.gz (Beagle reference)\n")
cat("- test_data/variants/genetic_map_chr11.txt (genetic map)\n")
cat("- test_data/samplesheet.csv (pipeline input)\n")
cat("- test_data/test.config (test configuration)\n")
cat("- test_data/setup_test_data.sh (setup script)\n")
cat("\nNext steps:\n")
cat("1. Run: cd test_data && ./setup_test_data.sh\n")
cat("2. Test pipeline: nextflow run main.nf -c test_data/test.config\n")
