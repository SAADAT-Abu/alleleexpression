#!/bin/bash -euo pipefail
bcftools view \
    --output P43_WT43.vcf.gz \
    -r chr11 \
    -f PASS -Oz \
    --threads 6 \
    P43_WT43.filtered.vcf.gz

tabix -p vcf P43_WT43.vcf.gz

cat <<-END_VERSIONS > versions.yml
"ASENEXT:BCFTOOLS_VIEW":
    bcftools: $(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*$//')
    tabix: $(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //; s/ .*$//')
END_VERSIONS
