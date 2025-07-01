#!/bin/bash -euo pipefail
printf "%s %s\n" 86500_ID2508_19-WT43T_S1_L001_R1_001.fastq.gz P43_WT43_1.gz 86500_ID2508_19-WT43T_S1_L001_R2_001.fastq.gz P43_WT43_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
    --quiet \
    --threads 6 \
    --memory 1365.3333333333 \
    P43_WT43_1.gz P43_WT43_2.gz

cat <<-END_VERSIONS > versions.yml
"ASENEXT:FASTQC":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
