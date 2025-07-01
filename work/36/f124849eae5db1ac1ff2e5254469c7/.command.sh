#!/bin/bash -euo pipefail
STAR \
    --genomeDir STAR_index \
    --readFilesIn input1/86500_ID2508_19-WT43T_S1_L001_R1_001.fastq.gz input2/86500_ID2508_19-WT43T_S1_L001_R2_001.fastq.gz \
    --readFilesCommand zcat \
    --runThreadN 8 \
    --outFileNamePrefix P43_WT43. \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile gencode.v47.primary_assembly.annotation.gtf \
    --outSAMattrRGline 'ID:P43_WT43'  'SM:P43_WT43' 'PL:illumina' \
    --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1 --waspOutputMode SAMtag --varVCFfile P43_WT43.GTonly.vcf \


if [ -f P43_WT43.Unmapped.out.mate1 ]; then
    mv P43_WT43.Unmapped.out.mate1 P43_WT43.unmapped_1.fastq
    gzip P43_WT43.unmapped_1.fastq
fi
if [ -f P43_WT43.Unmapped.out.mate2 ]; then
    mv P43_WT43.Unmapped.out.mate2 P43_WT43.unmapped_2.fastq
    gzip P43_WT43.unmapped_2.fastq
fi

cat <<-END_VERSIONS > versions.yml
"ASENEXT:STAR_ALIGN_WASP":
    star: $(STAR --version | sed -e "s/STAR_//g")
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
    gawk: $(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*$//')
END_VERSIONS
