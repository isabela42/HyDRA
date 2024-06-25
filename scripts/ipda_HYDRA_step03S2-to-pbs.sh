#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Jun 01, 2023
Last modified on Jun 23, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 03S2 (short-reads) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step03S2-to-pbs.sh -i "path/to/input/folder" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 1 -c 1 -w "00:30:00"

## Input:

-i <path/to/input/folder>   path/from/working/dir/to/HYDRA_step03S1_quality-check-after-correction_FastQC_DATE/
                            This folder should contain the HTML FastQC from step03S1 files.
-p <PBS stem>               Stem for PBS file names
-e <email>                  Email for PBS job
-m <INT>                    Memory INT required for PBS job in GB
-c <INT>                    Number of CPUS required for PBS job
-w <HH:MM:SS>               Clock walltime required for PBS job

## Output:

logfile.txt                 Logfile with commands executed and date
PBS files                   PBS file created

Pipeline description:

#   READ PRE-PROCESSING
#   ------------------------------------------------------------------------------------------------------------
#   01 [S] Quality check raw files (1FastQC, 2MultiQC)
#   01 [L] Quality check raw files (1NanoPlot, 2NanoComp, 3FastQC, 4MultiQC, 5Fasta_splitter)
#   02 [S] Short-reads correction (1Rcorrector) & Filter uncorrectable pair end reads (2FUPER)
#   02 [L] Trim reads of adapters (1Porechop), average-quality/headcrop (2Chopper) & poly(A) tails (3Cutadapt)
#-->03 [S] Quality check after correcting (1FastQC, 2MultiQC)
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
        i) input="${OPTARG}";;       # Input path to folder with FastQC HTML files
        p) pbs_stem="${OPTARG}";;    # Stem for PBS file names
        e) email="${OPTARG}";;       # Email for PBS job
        m) mem="${OPTARG}";;         # Memory required for PBS job
        c) ncpus="${OPTARG}";;       # Number of CPUS required for PBS job
        w) walltime="${OPTARG}";;    # Clock walltime required for PBS job
        h) Help ; exit;;             # Print Help and exit
        v) echo "${version}"; exit;; # Print version and exit
        ?) echo script usage: bash ipda_HYDRA_step03S2-to-pbs.sh -i path/to/input/folder -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step03S2-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

# NA

#................................................
#  Required modules, softwares and libraries
#................................................

## Load tools from HPC
# For more info, see
# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>

# None required

## Path to user-installed tools

# MultiQC:
# <https://multiqc.info/docs/getting_started/installation/>
multiqc="/working/lab_julietF/isabelaA/tools/anaconda3/bin/multiqc"

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step03S2_MultiQC="HYDRA_step03S2_quality-check-after-correction_MultiQC_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step03S2_MultiQC}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step03S2-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input folder of files:       ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${out_path_step03S2_MultiQC}"
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
echo "## Executing bash ipda_HYDRA_step03S2-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input folder of files:       ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${out_path_step03S2_MultiQC}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
echo "#!/bin/sh" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo "##########################################################################" >> ${pbs_stem}_${thislogdate}.pbs
echo "#" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Script:  ${pbs_stem}_${thislogdate}.pbs" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Version: v01" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Email:   ${email}" >> ${pbs_stem}_${thislogdate}.pbs
echo "#" >> ${pbs_stem}_${thislogdate}.pbs
echo "##########################################################################" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write PBS directives
echo "#PBS -N ${pbs_stem}_${thislogdate}" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -r n" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -m abe" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -M ${email}" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write directory setting
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Set main working directory" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo "## Change to main directory" >> ${pbs_stem}_${thislogdate}.pbs
echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write load modules
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "## Load tools from HPC"' >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "# No HPC modules required"' >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "## Path to user-installed tools"' >> ${pbs_stem}_${thislogdate}.pbs
echo "echo \"${multiqc}\"" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write PBS command lines
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Run step" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "## Run MultiQC at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "${multiqc} ${input} -o ${out_path_step03S2_MultiQC}" >> ${pbs_stem}_${thislogdate}.pbs

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 03S2 jobs) at
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
sed -i 's,${multiqc},'"${multiqc}"',g' "$logfile"
sed -i 's,${out_path_step03S2_MultiQC},'"${out_path_step03S2_MultiQC}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v