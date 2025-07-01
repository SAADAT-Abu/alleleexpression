#!/bin/bash -euo pipefail
beagle \
    gt=P43_WT43.vcf.gz \
    out=P43_WT43_beagle \
    ref=ALL.chr11.GRCh38.phased.chrprefix.vcf.gz \
    map=plink.chr11.GRCh38_chrprefix.map \
    nthreads=2

cat <<-END_VERSIONS > versions.yml
"ASENEXT:BEAGLE5_BEAGLE":
    beagle: $(echo $(beagle 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*$//'))
END_VERSIONS
