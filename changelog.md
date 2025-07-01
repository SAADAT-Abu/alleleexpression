# SAADAT-Abu/ASENext/: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.0.0](https://github.com/SAADAT-Abu/ASENext/)] - [30/06/2025]

Initial release of SAADAT-Abu/ASENext/, created with the [nf-core](https://nf-co.re/) template.

### Added

- Initial pipeline implementation
- STAR-WASP alignment for allele-aware mapping
- UMI-tools integration for molecular deduplication
- Beagle phasing for haplotype reconstruction
- phaser integration for ASE quantification
- FastQC for quality control
- MultiQC reporting
- Comprehensive documentation
- Singularity support
- Integration with nf-core infrastructure

### Pipeline features

- **Input**: Paired-end RNA-seq data with UMI barcodes and corresponding VCF files
- **Alignment**: STAR with WASP mode for unbiased allele-specific mapping
- **Filtering**: WASP filtering to remove reads with mapping bias
- **Deduplication**: UMI-based duplicate removal using UMI-tools
- **Phasing**: Beagle haplotype phasing with optional reference panels
- **ASE Analysis**: phaser for allele-specific expression quantification
- **QC**: Comprehensive quality control and reporting

### Supported platforms

- Local execution
- Container platforms (Singularity)

### Reference genomes

- Human (GRCh38, GRCh37)

### Added

- Core pipeline development
- Module integration and testing
- Documentation 
- Community feedback incorporation
