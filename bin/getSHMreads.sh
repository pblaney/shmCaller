#!/usr/bin/env bash
# Extracts mutated reads in region of BAM and performs local variant calling
# Author: pblaney

# Debugging settings
set -euo pipefail

# Set variables to hold user defined directory that contains all input BAM files and
# name of output merged BAM
normalBam=$1
tumorBam=$2
regionName=$3
region=$4
reference=$5
tumorSampleId=$(basename "${tumorBam}" | sed 's|\..*bam$||')

# get somatic SNVs per sample, these are all consensus non-germline SNVs assumed to be SHM (this can
# be confirmed with mutational signature analysis)
# get reads with SNVs at the sites, this eliminates reference identiy reads
# apply same read filter as IgCaller (exclude secondary[256], PCR/optical duplicates[1024])
# plus MAPQ > 20, stricter than IgCaller
# Then extract reads that overlap with SNVs and exclude reads with no mismatch to reference
echo "Pulling names of reads with SHM ..."
echo
samtools view -hb -F 256 -F 1024 -q 20 -e '[NM]!=0' "${tumorBam}" "${region}" \
	| java -jar /opt/jvarkit/dist/jvarkit.jar sam2tsv -R "${reference}" > "${tumorSampleId}.shm.sam2tsv.txt"

grep -v '^#' "${tumorSampleId}.shm.sam2tsv.txt" \
	| cut -f 1 \
	| uniq > "${tumorSampleId}.shm.read_qnames.txt"

# use list of read names to pull them from the bam and calculate the coverage for clone size estimation
samtools view -hb -N "${tumorSampleId}.shm.read_qnames.txt" "${tumorBam}" "${region}" > "${tumorSampleId}.shm.reads.bam"
samtools index "${tumorSampleId}.shm.reads.bam"
samtools coverage -r "${region}" -o "${tumorSampleId}.shm.reads.bam.depth.txt" "${tumorSampleId}.shm.reads.bam"

# Call SNVs and extract VAF information
echo
echo "Extracting SHM SNVs and SNPs ..."
echo 
bcftools mpileup --fasta-ref "${reference}" -a FORMAT/AD -Ov "${normalBam}" "${tumorSampleId}.shm.reads.bam" -r "${region}"  \
	| bcftools call -mA - \
	| bcftools norm -m- - \
	| bcftools view -i 'FORMAT/GT[1]="alt"' \
	| bcftools view -i 'FORMAT/GT[0]="alt"' \
	| bcftools +fill-tags - -Ov -o "${tumorSampleId}.${regionName}.snp.vcf" -- -t FORMAT/VAF

bcftools mpileup --fasta-ref "${reference}" -a FORMAT/AD -Ov "${normalBam}" "${tumorSampleId}.shm.reads.bam" -r "${region}" \
	| bcftools call -mA - \
	| bcftools view -i 'FORMAT/GT[1]="alt"' \
	| bcftools view -e 'FORMAT/GT[0]="alt"' \
	| bcftools norm -m- - \
	| bcftools +fill-tags - -Ov -o "${tumorSampleId}.${regionName}.shm.snv.vcf" -- -t FORMAT/VAF

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%AD\t%DP\t%VAF]\n' "${tumorSampleId}.${regionName}.shm.snv.vcf" \
	| cut -f 1-4,9-12> "${tumorSampleId}.${regionName}.shm.snv.txt"

echo
echo "Step complete ..."
echo 
