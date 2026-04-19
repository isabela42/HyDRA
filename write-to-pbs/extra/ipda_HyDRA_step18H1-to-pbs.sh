#!/bin/bash

version="1.0.1"
usage(){
echo "
Written by Isabela Almeida
Created on Dec 12, 2023
Last modified on May 15, 2025
Version: ${version}

Description: Write and submit PBS jobs for Step 18H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HyDRA_step18H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 10 -c 1 -w "400:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/hydra15H*_read-support_Samtools_DATE/folder/stem.fasta
                            of .spliced.fasta and .unspliced.fasta files (with full file name)  

                            Col2:
                            path/from/working/dir/to/hydra16H1_transcriptome-annotation_BedTools_DATE/stem/stem.bed
                            
                            Col3:
                            path/from/working/dir/to/protein-coding_reference-annotation.bed

                            It does not matter if same stem.fasta 
                            appears more than once on this input file.
                   
-p <PBS stem>               Stem for PBS file names
-e <email>                  Email for PBS job
-m <INT>                    Memory INT required for PBS job in GB
-c <INT>                    Number of CPUS required for PBS job
-w <HH:MM:SS>               Clock walltime required for PBS job

## Output:

logfile.txt                 Logfile with commands executed and date
PBS files                   PBS files created

Pipeline description:

#   READ PRE-PROCESSING
#   ------------------------------------------------------------------------------------------------------------
#   01 [S] Quality check raw files (1FastQC, 2MultiQC)
#   01 [L] Quality check raw files (1NanoPlot, 2NanoComp, 3FastQC, 4MultiQC, 5Fasta_splitter)
#   02 [S] Short-reads correction (1Rcorrector) & Filter uncorrectable pair end reads (2FUPER)
#   02 [L] Trim reads of adapters (1Porechop), average-quality/headcrop (2Chopper) & poly(A) tails (3Cutadapt)
#   03 [S] Quality check after correcting (1FastQC, 2MultiQC)
#   03 [L] Quality check after trimming (1NanoPlot, 2NanoComp, 3FastQC, 4MultiQC, 5Fasta_splitter)
#   04 [S] Trim reads of adapters (1Trimmomatic) & check ends/average quality (2BBDuk)
#   04 [L] Long-reads correction (1RopeBWT2, 2FMLRC2-convert, 3FMLRC2)
#   05 [S] Quality check after trimming (1FastQC, 2MultiQC)
#   05 [L] Read-lenght check (1Fasta_splitter)
#   06 [S] Filter rRNA reads (1BBDuk) & fix read headers (1in-house script)
#   06 [L] Filter rRNA reads (1BBDuk)
#   07 [S] Map reads to genome (1Bowtie2 to flag contamination, 1SAMtools)
#   07 [L] Map reads to genome (1Bowtie2 to flag contamination, 1SAMtools)
#   08 [S] Assess strandness (1RSeQC) & Assess replicates correlation (2DeepTools)
#   08 [L] Assess strandness (1RSeQC)
#   09 [S] Final quality check (1FastQC, 2MultiQC)
#   09 [L] Read-lenght check (1Fasta_splitter)
#   10 [H] Summary of quality control & read error correction steps

#   HYBRID DE NOVO TRANSCRIPTOME ASSEMBLY
#   ------------------------------------------------------------------------------------------------------------
#   11 [H] Hybrid de novo transcriptome assembly (1rnaSPAdes, 2Trinity)
#   12 [H] Raw assembly quality check (1BUSCO, 2TrinityStats, 3TransRate, 4TransRate)
#   13 [H] Transcriptome read representation (1GMAP & Bowtie2 for rnaSPAdes, 2Bowtie2 for Trinity assembly)
#   14 [H] Align transcripts to genome and assess splicing (1GMAP)
#   15 [H] Read support (1rnaSPAdes based, 2Trinity based - CDHit, Minimap2, Bowtie2, SAMtools)
#   16 [H] Annotation (1GMAP, 1sam2bed, 1BedTools)
#   17 [H] Summary of hybrid de novo transcriptome assembly steps

#   LONG NONCODING RNA DISCOVERY
#   ------------------------------------------------------------------------------------------------------------
#-->18 [H] Predict coding potential (1ezLncPred)
#   19 [H] Identify long noncoding RNAs (1FEELnc)
#   20 [H] Define lncRNAs (1Bash)
#   21 [H] Retrieve metrics, annotation and filter-out protein-coding overlaps (1BedTools, 2PBLAT)
#   22 [H] Summary of lncRNA discovery steps

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
"
}

## Display help message
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

#  ____       _   _   _                 
# / ___|  ___| |_| |_(_)_ __   __ _ ___ 
# \___ \ / _ \ __| __| | '_ \ / _` / __|
#  ___) |  __/ |_| |_| | | | | (_| \__ \
# |____/ \___|\__|\__|_|_| |_|\__, |___/
#                             |___/     

## Save execution ID
pid=`echo $$` #$BASHPID

## User name within your cluster environment
user=`whoami`

## Exit when any command fails
set -e

#................................................
#  Input parameters from command line
#................................................

## Get parameters from command line flags
while getopts "i:p:e:m:c:w:h:v" flag
do
    case "${flag}" in
        i) input="${OPTARG}";;       # Input files including path
        p) pbs_stem="${OPTARG}";;    # Stem for PBS file names
        e) email="${OPTARG}";;       # Email for PBS job
        m) mem="${OPTARG}";;         # Memory required for PBS job
        c) ncpus="${OPTARG}";;       # Number of CPUS required for PBS job
        w) walltime="${OPTARG}";;    # Clock walltime required for PBS job
        h) Help ; exit;;             # Print Help and exit
        v) echo "${version}"; exit;; # Print version and exit
        ?) echo script usage: bash ipda_HyDRA_step18H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
           exit;;
    esac
done

#................................................
#  Set Logfile stem
#................................................

## Set Logfile stem
# it contains all the executed commands with date/time;
# the output files general metrics (such as size);
# and memory/CPU usage for all executions
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
human_thislogdate=`date`
logfile=logfile_ipda_hydra18H1-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

# Protein-coding genes only, obtained using:
# WARNING! ONLY RUN IF NEEDED:
# grep 'transcript_biotype "protein_coding"' /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes.bed >> /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_protein-coding.bed
#genome_proteins=/working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_protein-coding.bed # hydra development
#genome_proteins=/working/lab_julietF/mainaB/ReferenceGenomes/ProteinCoding_TranscriptomeGRCh38rel79.bed # gems

### WARNING! # PLEK may have switched the labels "coding" and "non-coding".

#................................................
#  Required modules, softwares and libraries
#................................................

# GCC 7.5.0
module_gcc=gcc/7.5.0

# BedTools 2.29.0
module_bedtools=bedtools/2.29.0

# ezLncPred 1.4 - it doesn't work 
#module_ezLncPred=ezLncPred/1.4

# R 3.3.1
module_R=R/3.3.1

# Python3 - python virtual environment inside ezLncPred module
module_python=python/3.6.1

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step18H1_ezLncPred="hydra18H1_coding-potential_ezLncPred_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step18H1_ezLncPred}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HyDRA_step18H1-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input files:                 ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${out_path_step18H1_ezLncPred}"
echo "## logfile will be saved as:    ${logfile}"
echo

# ____  _             _   _                   
#/ ___|| |_ __ _ _ __| |_(_)_ __   __ _       
#\___ \| __/ _` | '__| __| | '_ \ / _` |      
# ___) | || (_| | |  | |_| | | | | (_| |_ _ _ 
#|____/ \__\__,_|_|   \__|_|_| |_|\__, (_|_|_)
#                                 |___/  

#................................................
#  Print Execution info to logfile
#................................................

exec &> "${logfile}"

date
echo "## Executing bash ipda_HyDRA_step18H1-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input files:                 ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${out_path_step18H1_ezLncPred}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#!/bin/sh" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Script:  ${pbs_stem}_${file}_${thislogdate}.pbs" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Version: v01" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Email:   ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#PBS -N ${pbs_stem}_${file}_${thislogdate}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#PBS -r n" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#PBS -m abe" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#PBS -M ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Set main working directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "## Change to main directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_gcc}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_R}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_python}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_bedtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "export OMP_NUM_THREADS=${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Retrieve high-confidence protein-coding transcripts (for comparison) at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; mkdir -p ${out_path_step18H1_ezLncPred}/${file} ; annotation=`grep ${path_file} ${input} | cut -f2` ; main_annotation=`echo "$(basename "${annotation}" | sed 's/\(.*\)\..*/\1/')"` ; genome_proteins=`grep "${path_file}" ${input} | cut -f3 | sort | uniq`; echo "bedtools intersect -wo -s -f 0.6 -r -a ${annotation} -b ${genome_proteins} >> ${out_path_step18H1_ezLncPred}/${file}/${main_annotation}.proteincoding.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run ezLncPred with CPC2 at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "ezLncPred -i ${path_file} -o ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2 CPC2" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run ezLncPred with CPAT at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "ezLncPred -i ${path_file} -o ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT CPAT -p Human" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run ezLncPred with CNCI at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "ezLncPred -i ${path_file} -o ${out_path_step18H1_ezLncPred}/${file}/${file}.CNCI CNCI -p ve" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run ezLncPred with PLEK at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
##cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "python ${plek} -fasta ${path_file} -out ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK PLEK -thread ${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
## How Maina did with PLEK model inside ezLncPred library: (above is from PLEK original publication)
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'n=`wc -l'" ${path_file} | cut -d' '"' -f1`; h=`echo "((${n}/2+1)/2)*2"|bc`; split -d -a1 -l${h}'" ${path_file} ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_; cat ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_0 ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_1 >> ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_A; cat ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_1 ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_0 >> ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_B; rm -f ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_0 ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_1" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "ezLncPred -i ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_A -o ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.PLEKa PLEK --thread ${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "ezLncPred -i ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_B -o ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.PLEKb PLEK --thread ${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "rm -f ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_A ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK_B" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Separate coding and noncoding transcripts [CPC2] at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f1,7,8 ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2 | grep -w "coding" | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2.proteincoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f1,7,8 ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2 | grep -w "noncoding" | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2.noncoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Separate coding and noncoding transcripts [CPAT] at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cat ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT | tr \"'\" '\t' | cut -f2,8 | awk 'NR>0{if(\$2>=0.5){print}}' | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT.proteincoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cat ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT | tr \"'\" '\t' | cut -f2,8 | awk 'NR>0{if(\$2<0.5){print}}' | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT.noncoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Separate coding and noncoding transcripts [CNCI] at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cat ${out_path_step18H1_ezLncPred}/${file}/${file}.CNCI/CNCI.index | tail -n +2 | cut -f1,2,3 | grep -Pe \"\tcoding\" | awk 'NR>0{printf \"%s\t%.2f\n\", \$1,\$3}' | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CNCI.proteincoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cat ${out_path_step18H1_ezLncPred}/${file}/${file}.CNCI/CNCI.index | tail -n +2 | cut -f1,2,3 | grep -Pe \"\tnoncoding\" | awk 'NR>0{printf \"%s\t%.2f\n\", \$1,\$3}' | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CNCI.noncoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Separate coding and noncoding transcripts [PLEK] at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f1,3 ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK* | tr -d '>' | grep \"^Coding\" | cut -f2- | sed 's/ len=/\tlen=/g' | sed 's/ path=/\tpath=/g' | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.proteincoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f1,3 ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK* | tr -d '>' | grep \"^Non-coding\" | cut -f2- | sed 's/ len=/\tlen=/g' | sed 's/ path=/\tpath=/g' | sort | uniq >> ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.noncoding" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Check protein coding transcripts count at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; annotation=`grep ${path_file} ${input} | cut -f2` ; echo "number=\`cut -f1 ${annotation} ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2.proteincoding | tail -n +2 | sort | uniq -c | grep -c \" 2 \"\` ; echo \"${out_path_step18H1_ezLncPred}/${file}/${file}\tCPC2 Annotated ptcoding classified as ptcoding\t${number}\" >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CPC2.proteincodingtest" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; annotation=`grep ${path_file} ${input} | cut -f2` ; echo "number=\`cut -f1 ${annotation} ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT.proteincoding | tail -n +2 | sort | uniq -c | grep -c \" 2 \"\` ; echo \"${out_path_step18H1_ezLncPred}/${file}/${file}\tCPAT Annotated ptcoding classified as ptcoding\t${number}\" >> ${out_path_step18H1_ezLncPred}/${file}/${file}.CPAT.proteincodingtest" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; annotation=`grep ${path_file} ${input} | cut -f2` ; echo "number=\`cut -f1 ${annotation} ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.proteincoding | tail -n +2 | sort | uniq -ci | grep -c \" 2 \"\` ; echo \"${out_path_step18H1_ezLncPred}/${file}/${file}\tPLEK Annotated ptcoding classified as ptcoding\t${number}\" >> ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.proteincodingtest" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; annotation=`grep ${path_file} ${input} | cut -f2` ; echo "number=\`cut -f1 ${annotation} ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.noncoding | tail -n +2 | sort | uniq -ci | grep -c \" 2 \"\` ; echo \"${out_path_step18H1_ezLncPred}/${file}/${file}\tPLEK Annotated ptcoding classified as noncoding\t${number}\" >> ${out_path_step18H1_ezLncPred}/${file}/${file}.PLEK.proteincodingtest" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
#ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n38 ${pbs} | head -n1)" ; echo "#                       $(tail -n35 ${pbs} | head -n1)" ; echo "#                       $(tail -n32 ${pbs} | head -n1)" ; echo "#                       $(tail -n29 ${pbs} | head -n1)" ; echo "#                       $(tail -n26 ${pbs} | head -n1)" ; echo "#                       $(tail -n25 ${pbs} | head -n1)" ; echo "#                       $(tail -n24 ${pbs} | head -n1)" ; echo "#                       $(tail -n23 ${pbs} | head -n1)" ; echo "#                       $(tail -n20 ${pbs} | head -n1)" ; echo "#                       $(tail -n19 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n15 ${pbs} | head -n1)" ; echo "#                       $(tail -n12 ${pbs} | head -n1)" ; echo "#                       $(tail -n11 ${pbs} | head -n1)" ; echo "#                       $(tail -n8 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n3 ${pbs} | head -n1)" ; echo "#                       $(tail -n2 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n26 ${pbs} | head -n1)" ; echo "#                       $(tail -n23 ${pbs} | head -n1)" ; echo "#                       $(tail -n20 ${pbs} | head -n1)" ; echo "#                       $(tail -n17 ${pbs} | head -n1)" ; echo "#                       $(tail -n14 ${pbs} | head -n1)" ; echo "#                       $(tail -n13 ${pbs} | head -n1)" ; echo "#                       $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n9 ${pbs} | head -n1)" ; echo "#                       $(tail -n6 ${pbs} | head -n1)" ; echo "#                       $(tail -n5 ${pbs} | head -n1)" ; echo "#                       $(tail -n2 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done


date ## Status of all user jobs (including HYDRA step 18H1 jobs) at
qstat -u "$user"

# This will remove $VARNAMES from output file with the actual $VARVALUE
# allowing for easily retracing commands
sed -i 's,${input},'"${input}"',g' "$logfile"
sed -i 's,${pbs_stem},'"${pbs_stem}"',g' "$logfile"
sed -i 's,${email},'"${email}"',g' "$logfile"
sed -i 's,${mem},'"${mem}"',g' "$logfile"
sed -i 's,${ncpus},'"${ncpus}"',g' "$logfile"
sed -i 's,${walltime},'"${walltime}"',g' "$logfile"
sed -i 's,${human_thislogdate},'"${human_thislogdate}"',g' "$logfile"
sed -i 's,${thislogdate},'"${thislogdate}"',g' "$logfile"
sed -i 's,${user},'"${user}"',g' "$logfile"
sed -i 's,${module_R},'"${module_R}"',g' "$logfile"
sed -i 's,${module_gcc},'"${module_gcc}"',g' "$logfile"
sed -i 's,${module_python},'"${module_python}"',g' "$logfile"
sed -i 's,${module_bedtools},'"${module_bedtools}"',g' "$logfile"
sed -i 's,${out_path_step18H1_ezLncPred},'"${out_path_step18H1_ezLncPred}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v
