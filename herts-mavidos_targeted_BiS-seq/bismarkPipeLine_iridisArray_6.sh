#!/bin/bash
# qsub parameters
#PBS -t 1-48 # job array by fastq file pair (n pairs) ##!!~~~ 
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=4
#PBS -l pmem=4gb
# File locations

# Bismark_v0.19.1 to v0.20.0 2018-09-05

# Usage instructions for script
usage() {
cat <<- EOF
bsmapPipieline_4.sh [options]

-h	help
-d	input directory
-r	runDir
-g	targeted genome
-G	whole genome
-w	writeover (doesn't do anything yet)
EOF
1>&2;
}

# Set woring directory
cd $PBS_O_WORKDIR
printf "## working directory: " 1>&2
pwd 1>&2

# Default values

INDIR="../data"
# # bsmap
# targetedGenome="/scratch/rja1e16/rat/genomes_Idx/targets/targetedCGIs.fa"
# wholeGenome="/scratch/rja1e16/rat/genomes_Idx/full/rn6.fa"
# # bismark
# targetedGenome="/scratch/rja1e16/rat/genomes_Idx/targets500bpFlank"
targetedGenome="/scratch/rja1e16/tRNA/targetedBSseq/genomes_Idx/"
wholeGenome="/scratch/rja1e16/genericResources/genomes/hg19/hg19-single/"

# target/WG/both?
# path to genome(s) -g and -G
# Argument parsing
# source location? - trim galore etc?

arrayRegex='([0-9]+).*'
[[ "$PBS_JOBID" =~ $arrayRegex ]]
runDir="run_${BASH_REMATCH[1]}"

while getopts hd:w:r:g:G: option; do
	case "${option}" in
		d) 
			INDIR=${OPTARG};;
		w) 
			WRITEOVER=${OPTARG};;
		r)
			runDir=${OPTARG};;
		g)
			targetedGenome=${OPTARG};;
		G)
			wholeGenome=${OPTARG};;
		h)
			usage
			exit 0;;
		*) 
			usage
			exit 1;;
	esac
done

# `date +%Y-%m-%d_%H-%M` run name? using job array id

if [[ ! -d "$runDir" ]]; then
	mkdir "$runDir"
fi
#echo ${BASH_REMATCH[1]}
#/home/rja1e16/software/trim_galore/
#/home/rja1e16/software/bismark/

#/rat
#- /src *script here
#- /data
#- /genomes_Idx
#- /fastqc-pre
#- /trimmed
#- /fastqc-post
#- /alignments
#- /methCalls

#export PERL5LIB="/home/rja1e16/perl5/lib/perl5/":$PERL5LIB


# Load modules from iridis
printf "## loading biobuilds/2017.11...\n" 1>&2
module load biobuilds/2017.11 #fastqc, samtools, bowtie2, cutadapt FASTQ
printf "## loaded\n" 1>&2
# Array of dirs containing paired fastqs
printf "## Getting data files...\n" 1>&2
#dirs=(../data/*)
dirs=("$INDIR"/*)
##!! dirs check?
printf "## working with: " 1>&2
run=${dirs[$PBS_ARRAYID]##*/}
printf "$run\n" 1>&2

# NB get fastq files  in dir
files=("$INDIR"/"$run"/*)

for i in "${!files[@]}"
do
   files[$i]="${files[$i]##*/}"
   files[$i]="${files[$i]%%.*}"
   printf "${files[$i]}\n"
done



########################################################################################
#### preliminary fastQC
########################################################################################
# fastqc-pre
printf "## making pre-trim fastqc dir...\n" 1>&2
if [[ ! -d "$runDir"/fastqc-pre ]]; then
	mkdir "$runDir"/fastqc-pre
fi
mkdir "$runDir"/fastqc-pre/"$run"

printf "## running pre-trim fastqc...\n" 1>&2


fastqc -o "$runDir"/fastqc-pre/"$run" \
"$INDIR"/"$run"/"${files[0]}"* \
"$INDIR"/"$run"/"${files[1]}"*

printf "## pre-trim fastqc completed\n" 1>&2

########################################################################################
#### Trimming
########################################################################################

############################################
## Fluidigm adaptor trimming 
############################################
# Adaptor and quality trimming
printf "## making trimmed reads dir...\n" 1>&2
if [[ ! -d "$runDir"/trimmed ]]; then
	mkdir "$runDir"/trimmed
fi
mkdir "$runDir"/trimmed/"$run"

printf "## trimming reads...\n" 1>&2

perl /scratch/rja1e16/tRNA/targetedBSseq/src/primerChecks_2_Linear.pl \
-p "$INDIR"/PrimerSetsCombined.txt \
-1 "$INDIR"/"$run"/"${files[0]}"* \
-2 "$INDIR"/"$run"/"${files[1]}"* \
-d "$runDir"/trimmed/"$run"/ \
--oCutadapt "$runDir"/trimmed/"$run"/"${files[0]}".fq \
--pCutadapt "$runDir"/trimmed/"$run"/"${files[1]}".fq \
--noseparate \
-t > "$runDir"/trimmed/"$run"/cutadaptTab.tab

printf "## first trimming completed...\n" 1>&2

# cutadapt \
#  -a 3p1RTc2=AGACCAAGTCTCTGCTACCGTA \
#  -A 3p2RTc1=TGTAGAACCATGTCGTCAGTGT \
#  -g 5p1c1=ACACTGACGACATGGTTCTACA \
#  -G 5p2c2=TACGGTAGCAGAGACTTGGTCT \
#  -o "$runDir"/trimmed/"$run"/"${files[0]}".fq \
#  -p "$runDir"/trimmed/"$run"/"${files[1]}".fq \
#  "$INDIR"/"$run"/"${files[0]}"* \
#  "$INDIR"/"$run"/"${files[1]}"* 

############################################
## Illumina Adaptor and quality trimming
############################################
printf "## making trimmed reads dir for custom primer trim...\n" 1>&2
if [[ ! -d "$runDir"/trimmed2 ]]; then
	mkdir "$runDir"/trimmed2
fi
mkdir "$runDir"/trimmed2/"$run"

/home/rja1e16/software/TrimGalore-0.5.0/trim_galore --paired \
 --trim1 \
 --max_n 1 \
 --quality 25 \
 -o "$runDir"/trimmed2/"$run" \
 "$runDir"/trimmed/"$run"/"${files[0]}"*.fq \
 "$runDir"/trimmed/"$run"/"${files[1]}"*.fq
printf "## second trimming completed...\n" 1>&2

########################################################################################
#### post-trim fastQC
########################################################################################

# # fastqc-post trimming
printf "## making post-trim fastqc dir...\n" 1>&2
if [[ ! -d "$runDir"/fastqc-post ]]; then
	mkdir "$runDir"/fastqc-post
fi
mkdir "$runDir"/fastqc-post/"$run"

printf "## running post-trim fastqc...\n" 1>&2
fastqc -o "$runDir"/fastqc-post/"$run" \
"$runDir"/trimmed2/"$run"/"${files[0]}"*.fq \
"$runDir"/trimmed2/"$run"/"${files[1]}"*.fq
printf "## post-trim fastqc completed\n" 1>&2

# ########################################################################################
# #### Alignments
# ########################################################################################

# ############################################
# ## Whole genome
# ############################################

# # Alignment
printf "## making alignment dir...\n" 1>&2

if [[ ! -d "$runDir"/alignments ]]; then
	mkdir "$runDir"/alignments
fi
mkdir "$runDir"/alignments/"$run"

printf "## aligning reads to whole genome...\n" 1>&2

# bsmap \
# -p 1 \
# -a "$runDir"/trimmed2/"$run"/"${files[0]}"*.fq \
# -b "$runDir"/trimmed2/"$run"/"${files[1]}"*.fq \
# -d "$wholeGenome" \
# -o "$runDir"/alignments/"$run"/"${files[0]}".bam \

/home/rja1e16/software/Bismark_v0.20.0/bismark -bowtie2 -n 1 \
 --genome "$wholeGenome" \
 --rg_tag \
 -o "$runDir"/alignments/"$run" \
 -1 "$runDir"/trimmed2/"$run"/"${files[0]}"*.fq \
 -2 "$runDir"/trimmed2/"$run"/"${files[1]}"*.fq 
printf "## alignment completed\n" 1>&2

# --non_directional \
############################################
## Targeted
############################################

# Targeted Alignment
##!! NB add flanks (~500bp?) to target regions !!##
printf "## making targeted alignment dir...\n" 1>&2
if [[ ! -d "$runDir"/targetAlignments ]]; then
	mkdir "$runDir"/targetAlignments
fi
mkdir "$runDir"/targetAlignments/"$run"
printf "## aligning reads target regions...\n" 1>&2

# bsmap \
# -p 1 \
# -a "$runDir"/trimmed2/"$run"/"${files[0]}"*.fq \
# -b "$runDir"/trimmed2/"$run"/"${files[1]}"*.fq \
# -d "$targetedGenome" \
# -o "$runDir"/alignments/"$run"/"${files[0]}".bam \

/home/rja1e16/software/Bismark_v0.20.0/bismark -bowtie2 -n 1 \
 --genome "$targetedGenome" \
 --rg_tag \
 -o "$runDir"/targetAlignments/"$run" \
 -1 "$runDir"/trimmed2/"$run"/"${files[0]}"*.fq \
 -2 "$runDir"/trimmed2/"$run"/"${files[1]}"*.fq 
# --rg_sample "$run" \ # needs --rg_id also
printf "## Target region alignment completed\n" 1>&2

########################################################################################
#### BSeQC - bias correction
########################################################################################

# printf "## BSeQC test...\n" 1>&2
# bseqc mbias \
# -s "$runDir"/alignments/"$run"/"${files[0]}"*.bam \
# -r /scratch/rja1e16/rat/genomes_Idx/full/rn6.fa \
# -n $run \
# -a
# printf "## BSeQC complete\n"

########################################################################################
#### Compressing fastq files
########################################################################################

# printf "## gzipping trimmed reads...\n" 1>&2
# compress trimmed alignments
# gzip "$runDir"/trimmed/"$run"/"${files[0]}"*.fq
# gzip "$runDir"/trimmed/"$run"/"${files[1]}"*.fq
# gzip "$runDir"/trimmed2/"$run"/"${files[0]}"*.fq
# gzip "$runDir"/trimmed2/"$run"/"${files[1]}"*.fq
# printf "## Trimmed reads compressed\n" 1>&2

########################################################################################
#### Methylation Calling
########################################################################################

############################################
## all
############################################
# Methylation calling - all
printf "## Making methCalls (all) dir...\n" 1>&2
if [[ ! -d "$runDir"/methCalls ]]; then
	mkdir "$runDir"/methCalls
fi
mkdir "$runDir"/methCalls/"$run"
# --scaffolds
printf "## Calling methylation levels...\n" 1>&2
/home/rja1e16/software/Bismark_v0.20.0/bismark_methylation_extractor \
 -p --no_overlap --gzip --scaffolds \
 --bedGraph \
 --genome_folder "$wholeGenome" \
 -o "$runDir"/methCalls/"$run"/ \
"$runDir"/alignments/"$run"/"${files[0]}"*.bam
printf "## Methylation calling completed\n" 1>&2
#--cytosine_report (2018-09-05) 
############################################
## Targeted
############################################

# Methylation calling - targeted
printf "## Making methCalls (targeted) dir...\n" 1>&2
if [[ ! -d "$runDir"/methCallsTargeted ]]; then
	mkdir "$runDir"/methCallsTargeted
fi
mkdir "$runDir"/methCallsTargeted/"$run"

printf "## Calling methylation levels...\n" 1>&2
/home/rja1e16/software/Bismark_v0.20.0/bismark_methylation_extractor \
 -p --no_overlap --gzip \
 --bedGraph \
 --genome_folder "$targetedGenome" \
 -o "$runDir"/methCallsTargeted/"$run"/ \
 "$runDir"/targetAlignments/"$run"/"${files[0]}"*.bam

# --cytosine_report
# --zero_based \
# --scaffolds \
printf "## Methylation calling completed\n" 1>&2

########################################################################################
#### reformatting bedgraph from targeted regions
########################################################################################
printf "## Reformating targeted bedgraph...\n"
# converting 'bedgraphs' from targeted alignment
cd $PBS_O_WORKDIR
bedgraph="$runDir"/methCallsTargeted/"$run"/*.bedGraph.gz
zcat $bedgraph | perl /scratch/rja1e16/rat/src/bedGraphTargetedExtractor.pl | gzip -c > "${bedgraph%/*}"/"$run".bedGraph.gz

coverage="$runDir"/methCallsTargeted/"$run"/*.bismark.cov.gz
zcat $coverage | perl /scratch/rja1e16/rat/src/bedGraphTargetedExtractor.pl | gzip -c > "${coverage%/*}"/"$run".bismark.cov.gz

printf "## Reformated targeted bedgraph\n"

