#!/usr/bin/env bash
# Extracts mutated reads in region of BAM and performs local variant calling
# Author: pblaney

# Debugging settings
set -uo pipefail

echo "######################################"
echo "#           Build SHM Clones         #"
echo "######################################"
echo

# Load modules
module add singularity/3.9.8
module add samtools/1.20
module add bcftools/1.20

# List modules for quick debugging
module list -t
echo 

# Set variable to hold input file with tumor/normal BAM files, region of interest name and genomic cooridnates
inputFile=$1
outputDir=$2

# Make the output directory if not already exists
mkdir -p "${outputDir}"

# Function to get SHM reads, call SHM SNVs / SNPs, cluster variants based on VAF density,
# extract the cluster specific reads, and build the SHM clones from the reads
buildShmClones() {
	while read -r line
		do
			# Collect input parameters from each line in input file
			normalBam=$(echo $line | cut -d ' ' -f 1)
			tumorBam=$(echo $line | cut -d ' ' -f 2)
			regionName=$(echo $line | cut -d ' ' -f 3)
			regionCoords=$(echo $line | cut -d ' ' -f 4)
			tumorSampleId=$(basename "${tumorBam}" | sed 's|\..*bam$||')

			# Header for separating samples runs in log
			echo "================================================================================"
			echo "SAMPLE: ${tumorSampleId}"
			echo 

			# Execute the getSHMreads script
			echo "Starting  >>>--->  getSHMreads"
			echo
			singularity exec -B $PWD:/temp/data --pwd /temp/data jvarkit.simg \
				./getSHMreads.sh \
				"normalBams/${normalBam}" \
				"tumorBams/${tumorBam}" \
				"${regionName}" \
				"${regionCoords}" \
				Homo_sapiens_assembly38.fasta

			sleep 2
			# Execute the cloneReadExtraction script
			echo "Starting  >>>--->  cloneReadExtraction"
			echo
			singularity exec -B $PWD:/temp/data --pwd /temp/data bmix-1.0.0.sif \
				./cloneReadExtraction.R \
				"${tumorSampleId}.${regionName}.shm.snv.txt" \
				"${tumorSampleId}.shm.sam2tsv.txt"

			sleep 2
			# Now build the clone specific BAMs and consensus FASTAs
			echo "Starting  >>>--->  buildClones"
			echo
			clonesDetected=$(ls -1 ${tumorSampleId}_Bin*_read_qnames.txt | wc -l)
			echo "Found ${clonesDetected} clones"
			for reads in `ls -1 ${tumorSampleId}_Bin*_read_qnames.txt`
			do
				# Create names of output files
				clusterId=$(echo ${reads} | sed 's|_read_qnames\.txt$||')
				echo "* ${clusterId}"

				# Build the per-cluster BAMs using file of names and calculate the coverage for clone size estimation
				samtools view -hb -N "${reads}" "tumorBams/${tumorBam}" "${regionCoords}" > "${clusterId}.bam"
				samtools index "${clusterId}.bam"
				samtools coverage -r "${regionCoords}" -o "${clusterId}.bam.depth.txt" "${clusterId}.bam"

				samtools consensus -f fasta -a --show-del yes --show-ins yes "${clusterId}.bam" -r "${regionCoords}" \
					| sed "s|^>.*|>${clusterId}_clone|" > "${clusterId}.fasta"
			done

			# Next, merge all clone consensus FASTAs in to single file
			echo "Merging all clone-specific FASTAs to consensus FASTA for MSA analysis"
			cat ${tumorSampleId}*.fasta > "${tumorSampleId}.consensus.fasta"
			
			# Finally build the BAM of the reads not included in the complete set of per-cluster BAMs
			echo "Collecting for null bin reads"
			grep -Fxv -f <(cat ${tumorSampleId}_Bin*_read_qnames.txt) ${tumorSampleId}.shm.read_qnames.txt > "${tumorSampleId}_nullBin_read_qnames.txt"

			if [[ -s "${tumorSampleId}_nullBin_read_qnames.txt" ]]; then
				echo "* Found nullBin reads"
				samtools view -hb -N "${tumorSampleId}_nullBin_read_qnames.txt" "tumorBams/${tumorBam}" "${regionCoords}" > "${tumorSampleId}_nullBin.bam"
				samtools index "${tumorSampleId}_nullBin.bam"
				samtools coverage -r "${regionCoords}" -o "${tumorSampleId}_nullBin.bam.depth.txt" "${tumorSampleId}_nullBin.bam"
			else
				echo "WARNING: Zero nullBin reads found, confirm this by checking clones depth with SHM BAM depth"
			fi

			echo
			echo "Step complete ..."
			echo 

			sleep 2
			# Move all output files to output directory
			mkdir -p "${outputDir}${tumorSampleId}/"
			mv ${tumorSampleId}* "${outputDir}${tumorSampleId}/"

			# Footer for each run
			echo "${tumorSampleId} ... DONE"

		done < "${inputFile}"
}

# Call the function
buildShmClones

echo "================================================================================"
echo
