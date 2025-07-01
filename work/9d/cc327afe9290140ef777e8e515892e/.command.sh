#!/bin/bash -euo pipefail
# Create index if it doesn't exist
if [ ! -f "WT43.deepvariant.vcf.gz.tbi" ]; then
tabix -p vcf "WT43.deepvariant.vcf.gz"
fi
# Create GT-only uncompressed VCF for STAR
bcftools view -f PASS "WT43.deepvariant.vcf.gz" | bcftools annotate -x INFO,^FORMAT/GT -Ov -o "P43_WT43.GTonly.vcf"

# Create PASS-filtered compressed and indexed VCF for Beagle
bcftools view -r chr11 -f PASS -Oz -o "P43_WT43.filtered.vcf.gz" "WT43.deepvariant.vcf.gz"
tabix -p vcf "P43_WT43.filtered.vcf.gz"

cat <<-END_VERSIONS > versions.yml
"ASENEXT:PREPARE_VCF":
    bcftools: $(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*$//')
    tabix: $(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*$//')
END_VERSIONS
