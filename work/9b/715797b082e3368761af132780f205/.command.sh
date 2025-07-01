#!/bin/bash -euo pipefail
check_samplesheet.py \
    samplesheet.csv \
    samplesheet.valid.csv

cat <<-END_VERSIONS > versions.yml
"ASENEXT:SAMPLESHEET_CHECK":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
