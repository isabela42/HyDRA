#!/bin/bash

version="1.0.1"
usage(){
echo "
Written by Isabela Almeida
Created on Sep 12, 2023
Last modified on May 15, 2025
Version: ${version}

Description: Write and submit PBS jobs for Step 11H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HyDRA_step11H1-to-pbs.sh -i "path/to/input/folder" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 250 -c 40 -w "60:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                        
                            Col1:
                            path/from/working/dir/to/hydra06L1_ribodepletion_BBDuk_DATE/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed*
                            of ALL .fasta files and no full stops.
                            Please provide the stem that includes all files.

                            Col2:
                            path/from/working/dir/to/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm*
                            of ALL .cor.fq.gz files (corresponding to _R1 and _R2 files) and no full stops.
                            If different samples are included, please provide the stem that includes all of them

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
#-->11 [H] Hybrid de novo transcriptome assembly (1rnaSPAdes, 2Trinity)
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
        ?) echo script usage: bash ipda_HyDRA_step11H1-to-pbs.sh -i path/to/input/folder -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_hydra11H1-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# rnaSPAdes 3.14.1:
# <https://cab.spbu.ru/files/release3.14.1/rnaspades_manual.html>
# /software/SPAdes/SPAdes-3.14.1/bin/rnaspades.py 
module_rnaSPAdes="SPAdes/3.14.1"

# seqtk 1.3:
module_seqtk=seqtk/20191028

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step11H1_rnaSPAdes="hydra11H1_assembly_rnaSPAdes_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step11H1_rnaSPAdes}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HyDRA_step11H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step11H1_rnaSPAdes}"
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
echo "## Executing bash ipda_HyDRA_step11H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step11H1_rnaSPAdes}"
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
echo "module load ${module_rnaSPAdes}" >> ${pbs_stem}_${thislogdate}.pbs
echo "module load ${module_seqtk}" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write PBS command lines
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Run step" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "## Run gzip at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
pairedend=`cut -f2 ${input}`; longreads=`cut -f1 ${input}`; fileshort=`echo "$(basename "${pairedend%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; filelong=`echo "$(basename "${longreads%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`
echo "zcat -c ${pairedend}*_R1.cor.fq.gz | gzip -c >> ${out_path_step11H1_rnaSPAdes}/${fileshort}_merged-paired-end_R1.cor.fq.gz" >> ${pbs_stem}_${thislogdate}.pbs
echo "zcat -c ${pairedend}*_R2.cor.fq.gz | gzip -c >> ${out_path_step11H1_rnaSPAdes}/${fileshort}_merged-paired-end_R2.cor.fq.gz" >> ${pbs_stem}_${thislogdate}.pbs
echo "gzip -c ${longreads}*.fasta >> ${out_path_step11H1_rnaSPAdes}/${filelong}_merged.fa.gz" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Run rnaSPAdes at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "rnaspades.py --checkpoint "all" -t ${ncpus} -m ${mem} -1 ${out_path_step11H1_rnaSPAdes}/${fileshort}_merged-paired-end_R1.cor.fq.gz -2 ${out_path_step11H1_rnaSPAdes}/${fileshort}_merged-paired-end_R2.cor.fq.gz --nanopore ${out_path_step11H1_rnaSPAdes}/${filelong}_merged.fa.gz -o ${out_path_step11H1_rnaSPAdes}" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Get singleline fasta at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "seqtk seq ${out_path_step11H1_rnaSPAdes}/transcripts.fasta > ${out_path_step11H1_rnaSPAdes}/transcripts.singleline.fasta " >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Remove temporary files at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "rm -f ${out_path_step11H1_rnaSPAdes}/${fileshort}_merged-paired-end_R1.cor.fq.gz ${out_path_step11H1_rnaSPAdes}/${fileshort}_merged-paired-end_R2.cor.fq.gz ${out_path_step11H1_rnaSPAdes}/${filelong}_merged.fa.gz" >> ${pbs_stem}_${thislogdate}.pbs

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n12 ${pbs} | head -n1)" ; echo "                        $(tail -n11 ${pbs} | head -n1)" ; echo "                        $(tail -n10 ${pbs} | head -n1)" ; echo "                        $(tail -n7 ${pbs} | head -n1)" ; echo "                        $(tail -n4 ${pbs} | head -n1)" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 11H1 jobs) at
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
sed -i 's,${module_rnaSPAdes},'"${module_rnaSPAdes}"',g' "$logfile" 
sed -i 's,${module_seqtk},'"${module_seqtk}"',g' "$logfile"
sed -i 's,${out_path_step11H1_rnaSPAdes},'"${out_path_step11H1_rnaSPAdes}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v