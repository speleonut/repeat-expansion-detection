#!/bin/bash

#SBATCH -J ExpHunter
#SBATCH -o /hpcfs/users/%u/log/ExpHunt.slurm-%j.out

#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=01:00:00
#SBATCH --mem=8GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mark.corbett@adelaide.edu.au

neuroDir=/hpcfs/groups/phoenix-hpc-neurogenetics
refDir=$neuroDir/RefSeq
userDir=/hpcfs/users/${USER}
expHunterPath=$neuroDir/executables/ExpansionHunter_latest/bin

modList=("SAMtools/1.17-GCC-11.2.0")

usage()
{
echo "# Script for finding repeat expansions in Illumina short-read GS
# Requires: Text file with a list of bam or cram files including the full path to each file. Expansion Hunter v5+.
#
# Usage sbatch --array 0-(n-1 bam files) $0 -i inputFile.txt [-o /path/to/output] | [ - h | --help ]
#
# Options
# -i <arg>    REQUIRED. Text file with a list of bam or cram files including the full path to each file.
# -g <arg>    OPTIONAL. The script will try to locate the right genome based on the @SQ lines in the bam header if you don't set this. 
#                       Letting the script do the work is probably best. With this option you can specify the path to the original reference that your BAM or CRAM file was mapped to.
# -c <arg>    OPTIONAL. The script will select a full variant catalog matched to your genome. You can use this option to use a custom catalog.
# -o <arg>    OPTIONAL. Path to where you want to find your file output (if not specified $userDir/expansionHunter/output/prefix is used).
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
# Original: Mark Corbett, 12/02/2018
# Modified: (Date; Name; Description)
# 21/09/2018; Mark Corbett; Changed script to accept batch input
# 25/06/2020; Mark Corbett; Update to EH v3+
# 24/11/2020; Mark Corbett; Update the EH v4+
# 15/05/2024; Mark Corbett; Update to EH v5+, automate catalog selection
#
"
}

select_genome_build()
{
case "${genomeSize}" in
    3099922541 )    buildID="GRCh38"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
                    ;;
    3217346917 )    buildID="hs38DH"
                    genomeBuild="$refDir/hs38DH.fa"
                    ;;
    3137454505 )    buildID="hs37d5"
                    genomeBuild="$refDir/hs37d5.fa.gz"
                    ;;
    2730871774 )    buildID="GRCm38"
                    genomeBuild="$refDir/GRCm38_68.fa"
                    ;;
    3117463893 )    buildID="CHM13v2"
                    genomeBuild="$refDir/T2T_CHM13v2.0.ucsc.ebv.fa.gz"
                    ;;
    3137161264 )    buildID="hg19"
                    genomeBuild="$refDir/ucsc.hg19.fasta"
                    ;;
    3105715063 )    buildID="GRCh38.hs38d1"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
                    ;;
    3099750718 )    buildID="GRCh38"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
                    ;;
    3031042417 )    buildID="GRCh38.blacklist"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
                    ;;
    3101804741 )    buildID="hg19_1stM_unmask_ran_all"
                    genomeBuild="$refDir/hg19_1stM_unmask_ran_all.fa"
                    ;;
    * )         echo "## ERROR: Genome length $genomeSize for ${bamFile[SLURM_ARRAY_TASK_ID]} was not matched, you may need to specify the genome build directly using the -g flag."
                exit 1
                ;;
esac
echo "## INFO: The file ${bamFile[SLURM_ARRAY_TASK_ID]} was likely mapped to $buildID corresponding to the refseq $genomeBuild."

if [ -z "$repeatSpecs" ]; then # If the variant catalog was not specified then find one that matches the genome build
    case $buildID in
        GRCh38 | hs38DH | GRCh38.hs38d1 | GRCh38.blacklist )    repeatSpecs="$neuroDir/executables/ExpansionHunter_gnomAD_catalogs/variant_catalog_with_offtargets.GRCh38.json"
                                                                ;;
        hg19 | hg19_1stM_unmask_ran_all | hs37d5 )              repeatSpecs="$neuroDir/executables/ExpansionHunter_gnomAD_catalogs/variant_catalog_with_offtargets.GRCh37.json"
                                                                ;;
        CHM13v2 )                                               repeatSpecs="$neuroDir/executables/ExpansionHunter_gnomAD_catalogs/variant_catalog_without_offtargets.CHM13v2.json"
                                                                ;;
    esac
    echo "## INFO: Using the following variant catalog: $repeatSpecs"
fi 
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -i )                    shift
                                        inputFile=$1
                                        ;;
                -o )                    shift
                                        outDir=$1
                                        ;;
                -g )                    shift
                                        genomeBuild=$1
                                        ;;
                -c )                    shift
                                        repeatSpecs=$1
                                        ;;
                -h | --help )           usage
                                        exit 0
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

if [ -z "$inputFile" ]; then # If no input file specified then do not proceed
        usage
        echo "#ERROR: You need to give me a text file with a list of BAM or CRAM files including the full path to each file."
        exit 1
fi

readarray -t bamFile < $inputFile

# Load modules
for mod in "${modList[@]}"; do
    module load $mod
done

sampleID = $( samtools samples ${bamFile[SLURM_ARRAY_TASK_ID]} | cut -f1 )

if [ -z "$genomeBuild" ]; then # If genome not specified then see if it is possible to find the reference.  This will also match a variant catalog if a custom one was not specified.
    genomeSize=$(samtools view -H ${bamFile[SLURM_ARRAY_TASK_ID]} | grep @SQ | cut -f3 | cut -f2 -d":" | awk '{s+=$1} END {printf "%.0f\n", s}' -)
    select_genome_build
fi

if [ -z "$outDir" ]; then # If no output directory then use a default directory
        outDir=$userDir/expansionHunter/output/$sampleID
        echo "## INFO: Using $outDir as the output directory"
fi

# Make sure $outDir exists
if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

# Start the script
$expHunterPath/ExpansionHunter \
--reads ${bamFile[$SLURM_ARRAY_TASK_ID]} \
--reference $genomeBuild \
--variant-catalog $repeatSpecs \
--output-prefix $outDir/$sampleID
