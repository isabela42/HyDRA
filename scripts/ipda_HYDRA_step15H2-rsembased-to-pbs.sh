#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Oct 22, 2023
Last modified on Jun 23, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 15H2 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step15H2-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 80 -c 16 -w "320:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/HYDRA_step14H1_splicing-assessment_GMAP_DATE/stem.fasta
                            of .spliced.fasta and .unspliced.fasta files (with full file name)

                            Col2:
                            Splicing info, i.e.:
                            for path/to/file.spliced.fasta, Col2 should contain spliced
                            for path/to/file.unspliced.fasta, Col2 should contain unspliced

                            Col3:
                            path/from/working/dir/to/HYDRA_step11H2_de-novo-transcriptome-assembly_*_DATE/insilico_read_normalization/left.norm.fq.paired.fq
                            of .Trinity.fasta derived files.

                            Col4:
                            path/from/working/dir/to/HYDRA_step11H2_de-novo-transcriptome-assembly_*_DATE/insilico_read_normalization/right.norm.fq.paired.fq
                            of .Trinity.fasta derived files.

                            Col5:
                            INT of high FPKM cutoff for that file, e.g.: 1 for spliced and 5 for unspliced

                            Col6:
                            INT of low FPKM cutoff for that file, e.g.: 0.5 for spliced and 3 for unspliced                          

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
#-->15 [H] Read support (1rnaSPAdes based, 2Trinity based - CDHit, Minimap2, Bowtie2, SAMtools)
#   16 [H] Annotation (1GMAP, 1sam2bed, 1BedTools)
#   17 [H] Summary of hybrid de novo transcriptome assembly steps

#   LONG NONCODING RNA DISCOVERY
#   ------------------------------------------------------------------------------------------------------------
#   18 [H] Predict coding potential (1ezLncPred)
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
        ?) echo script usage: bash ipda_HYDRA_step15H2-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step15H2-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

# Memory required (in MB):
memreqmb=`echo "${mem}*1000" | bc`

# Memory per thread:
mempt=`echo "${mem}/${ncpus}" | bc`

# FPKM recommended filter cutoffs:
# high_fpkm_cut_spliced=1; high_fpkm_cut_unspliced=5
# low_fpkm_cut_spliced=0.5; low_fpkm_cut_unspliced=3

# CDhit cutoffs:
identity=98; coverage=90

#................................................
#  Required modules, softwares and libraries
#................................................

## Load tools from HPC
# For more info, see
# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>

# Trinity
# <https://github.com/trinityrnaseq/>
module_trinityold=trinityrnaseq/2.2.0
trinfilterold="/software/trinityrnaseq/trinityrnaseq-2.2.0/util/filter_fasta_by_rsem_values.pl"

# Bowtie2
module_bowtie2=bowtie2/2.2.9

# SamTools
module_samtools=samtools/1.9

# RSEM 1.3.1
module_rsem=RSEM/1.3.1

# cdhit 4.7
module_cdhit=cd-hit/4.6.8-2017-1208

# R 4.3.1
module_R=R/4.3.1

## Path to user-installed tools

# None required

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step15H2_RSEM="HYDRA_step15H2_FPKM-filter_RSEM_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step15H2_RSEM}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step15H2-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step15H2_RSEM}"
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
echo "## Executing bash ipda_HYDRA_step15H2-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step15H2_RSEM}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#!/bin/sh" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Script:  ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Version: v01" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Email:   ${email}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#PBS -N ${pbs_stem}_${file}_${splice_info}_${thislogdate}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#PBS -r n" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#PBS -m abe" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#PBS -M ${email}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#................................................" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Set main working directory" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#................................................" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "## Change to main directory" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#................................................" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#................................................" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Load tools from HPC"' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "module load ${module_trinityold}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "module load ${module_bowtie2}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "module load ${module_samtools}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "module load ${module_rsem}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "module load ${module_cdhit}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "module load ${module_R}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Path to user-installed tools"' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "# None required"' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#................................................" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#  Run step" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "#................................................" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Filter with CDhit based on sequence similarity at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; mkdir -p  ${out_path_step15H2_RSEM}/${file}_${splice_info} ; echo "cd-hit-est -o ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_CDhit_identity${identity}_cov${coverage}.${splice_info}.fasta -c 0.${identity} -i ${path_file} -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T ${ncpus} -M ${memreqmb}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Build transcriptome bowtie2 index at at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "bowtie2-build --threads ${ncpus} ${path_file} ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_bt2idx" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Run Bowtie2 alignment at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; f1=`grep ${path_file} ${input} | cut -f3` ; f2=`grep ${path_file} ${input} | cut -f4`; echo "bowtie2 -q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --threads ${ncpus} -x ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_bt2idx -1 ${f1} -2 ${f2} | samtools view -@${ncpus} -bS - >> ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.bam" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Coordination sort transcriptome alignment at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "samtools sort -@${ncpus} ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.bam -o ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.coordsorted.bam" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "rm -rf ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.bam" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Run RSEM prepare transcriptome idx at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "rsem-prepare-reference --num-threads ${ncpus} ${path_file} ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_rsemidx" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Run convert-sam-for-rsem at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "convert-sam-for-rsem --num-threads ${ncpus} ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.coordsorted.bam ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.coordsorted.rsem" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Run RSEM calculate expression at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "rsem-calculate-expression --num-threads ${ncpus} --no-bam-output --alignments --paired-end ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads_bt2.coordsorted.rsem.bam ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_rsemidx ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads-rsem-expr" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Filter high FPKM based on RSEM results at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; high_fpkm_cutoff=`grep ${path_file} ${input} | cut -f5` ; echo "${trinfilterold} --rsem_output=${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads-rsem-expr.isoforms.results --fasta ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_CDhit_identity${identity}_cov${coverage}.${splice_info}.fasta --output ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_shortreads-rsem-expr_FPKM_highcutoff${high_fpkm_cutoff}.${splice_info}.fasta --fpkm_cutoff=${high_fpkm_cutoff}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Filter low FPKM based on RSEM results at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; low_fpkm_cutoff=`grep ${path_file} ${input} | cut -f6` ; echo "${trinfilterold} --rsem_output=${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads-rsem-expr.isoforms.results --fasta ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_CDhit_identity${identity}_cov${coverage}.${splice_info}.fasta --output ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_shortreads-rsem-expr_FPKM_lowcutoff${low_fpkm_cutoff}.${splice_info}.fasta --fpkm_cutoff=${low_fpkm_cutoff}" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo 'echo "## Run RSEM plot-model at" ; date ; echo' >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; splice_info=`grep ${path_file} ${input} | cut -f2` ; echo "rsem-plot-model ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads-rsem-expr ${out_path_step15H2_RSEM}/${file}_${splice_info}/${file}_${splice_info}_shortreads-rsem-expr.model.pdf" >> ${pbs_stem}_${file}_${splice_info}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n29 ${pbs} | head -n1)" ; echo "#                       $(tail -n26 ${pbs} | head -n1)" ; echo "#                       $(tail -n23 ${pbs} | head -n1)" ; echo "#                       $(tail -n20 ${pbs} | head -n1)" ; echo "#                       $(tail -n19 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n13 ${pbs} | head -n1)" ; echo "#                       $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 15H2 jobs) at
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
sed -i 's,${module_trinityold},'"${module_trinityold}"',g' "$logfile"
sed -i 's,${module_bowtie2},'"${module_bowtie2}"',g' "$logfile"
sed -i 's,${module_samtools},'"${module_samtools}"',g' "$logfile"
sed -i 's,${module_rsem},'"${module_rsem}"',g' "$logfile"
sed -i 's,${module_cdhit},'"${module_cdhit}"',g' "$logfile"
sed -i 's,${module_R},'"${module_R}"',g' "$logfile" 
sed -i 's,${memreqmb},'"${memreqmb}"',g' "$logfile"
sed -i 's,${mempt},'"${mempt}"',g' "$logfile"
sed -i 's,${identity},'"${identity}"',g' "$logfile"
sed -i 's,${coverage},'"${coverage}"',g' "$logfile"
sed -i 's,${out_path_step15H2_RSEM},'"${out_path_step15H2_RSEM}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v
