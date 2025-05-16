#!/bin/bash

version="1.0.1"
usage(){
echo "
Written by Isabela Almeida
Created on Feb 14, 2024
Last modified on May 15, 2025
Version: ${version}

Description: Write and submit PBS jobs for Step 21H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HyDRA_step21H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 1 -c 3 -w "01:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/hydra18H1_coding-potential_ezLncPred_DATE/stem/
                            ezlncpred-and-feelnc.noncoding

                            Col2:
                            path/from/working/dir/to/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/
                            stem.gtf

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
#   18 [H] Predict coding potential (1ezLncPred)
#   19 [H] Identify long noncoding RNAs (1FEELnc)
#   20 [H] Define lncRNAs (1Bash)
#-->21 [H] Retrieve metrics, annotation and filter-out protein-coding overlaps (1BedTools, 2PBLAT)
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
        ?) echo script usage: bash ipda_HyDRA_step21H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_hydra21H1-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# Perl 5.22
module_perl=perl/5.22

# HTSlib 1.19.1
module_htslib=htslib/1.19.1

# Bedtools 2.29.0
module_bedtools=bedtools/2.29.0 

# gtf2gff 0.1
# <https://manpages.ubuntu.com/manpages/trusty/man1/gtf2gff3.1p.html>
gtf2gff=/working/lab_julietF/isabelaA/tools/gtf2gff3.pl

# genestats 1.0
genestats=/working/lab_julietF/mainaB/Tools/genestats.sh

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step21H1_Bash="hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step21H1_Bash}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HyDRA_step21H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step21H1_Bash}"
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
echo "## Executing bash ipda_HyDRA_step21H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step21H1_Bash}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#!/bin/sh" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Script:  ${pbs_stem}_${file}_${thislogdate}.pbs" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Version: v01" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Email:   ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#PBS -N ${pbs_stem}_${file}_${thislogdate}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#PBS -r n" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#PBS -m abe" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#PBS -M ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Set main working directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "## Change to main directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "module load ${module_perl}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "module load ${module_htslib}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "module load ${module_bedtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#${gtf2gff}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#${genestats}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo 'echo "## Get lncRNA GTF at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; mkdir -p ${out_path_step21H1_Bash}/${file} ; minimap=`grep ${path_file} ${input} | cut -f2` ; echo "cat ${path_file}"' | while read transcript; do grep "gene_id \"${transcript}\";"'" ${minimap} >> ${out_path_step21H1_Bash}/${file}/${file}.gtf; done" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo 'echo "## Get lncRNA GFF at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "perl ${gtf2gff} ${out_path_step21H1_Bash}/${file}/${file}.gtf >> ${out_path_step21H1_Bash}/${file}/${file}.gff" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo 'echo "## Get lncRNA genestats at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "${genestats} ${out_path_step21H1_Bash}/${file}/${file}.gff >> ${out_path_step21H1_Bash}/${file}/${file}.genestats" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo 'echo "## Get lncRNA metrics at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=$(echo "$path_file" | awk -F'/' '{print $NF}' | awk -F'.' '{for (i=1; i<NF-2; i++) printf "%s.", $i; print $(NF-2)}') ; echo "cat ${out_path_step21H1_Bash}/${file}/${file}.genestats | while read line; do echo"' "${line}" | awk -v OFS="\t" '"'"'{ exon_sum += $3; exon_len += $4; intron_sum += $5; intron_len += $6 } END { print exon_sum, exon_len / exon_sum, intron_sum, intron_len / (intron_sum+1) }'"' >> ${out_path_step21H1_Bash}/${file}/${file}.metrics ; done" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 21H1 jobs) at
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
sed -i 's,${module_perl},'"${module_perl}"',g' "$logfile"
sed -i 's,${module_htslib},'"${module_htslib}"',g' "$logfile"
sed -i 's,${module_bedtools},'"${module_bedtools}"',g' "$logfile"
sed -i 's,${minimap},'"${minimap}"',g' "$logfile" 
sed -i 's,${gtf2gff},'"${gtf2gff}"',g' "$logfile"
sed -i 's,${genestats},'"${genestats}"',g' "$logfile"
sed -i 's,${out_path_step21H1_Bash},'"${out_path_step21H1_Bash}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v