#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Feb 13, 2024
Last modified on Jun 23, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 22H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step22H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 1 -c 1 -w "00:30:00"

## Input:

-i <path/to/input/files>    TXT file containing on:
                            Col1: Folder run name (no path)
                            Col2: Path/to/logfile
                            Col3: Path/to/pbs/*_DATE.e*
                            Col4: Path/to/pbs/*_DATE.o*
                            Col5: Path/to/pbs/*_DATE.pbs
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
#   21 [H] Retrieve metrics, annotation and filter-out protein-coding overlaps (1BedTools, 2PBLAT)
#-->22 [H] Summary of lncRNA discovery steps

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
        ?) echo script usage: bash ipda_HYDRA_step22H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step22H1-to-pbs_${thislogdate}.txt

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

# None required

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step22H1_Bash="HYDRA_step22H1_summary_Bash_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step22H1_Bash}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step22H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step22H1_Bash}"
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
echo "## Executing bash ipda_HYDRA_step22H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step22H1_Bash}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Summary 18H1
#................................................

date ## Collect run summary for step 18H1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step18H1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tInput\tCNCI-noncoding\tCNCI-coding\tCPAT-noncoding\tCPAT-coding\tCPC2-noncoding\tCPC2-coding\tMin2\tAll\tBase-file"
summary_output="${out_path_step22H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

fasta=`head -n1 ${pbs_e} | cut -d" " -f4 | tr -d "'"`;
main=`echo "$(basename ${fasta})" | sed 's/.\{6\}$//'`;

col01="${step}";                                                           # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f5`;                       # Sample
col03=`grep -c "^>" ${fasta}`;                                             # Input 
col04=`awk 'END{print NR}' ${folder}/${main}/${main}.CNCI.noncoding`;      # CNCI-noncoding
col05=`awk 'END{print NR}' ${folder}/${main}/${main}.CNCI.proteincoding`;  # CNCI-coding
col06=`awk 'END{print NR}' ${folder}/${main}/${main}.CPAT.noncoding`;      # CPAT-noncoding
col07=`awk 'END{print NR}' ${folder}/${main}/${main}.CPAT.proteincoding`;  # CPAT-coding
col08=`awk 'END{print NR}' ${folder}/${main}/${main}.CPC2.noncoding`;      # CPC2-noncoding
col09=`awk 'END{print NR}' ${folder}/${main}/${main}.CPC2.proteincoding`;  # CPC2-coding
col10=`echo "$(basename ${pbs_e})"`;                                       # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}" >> ${summary_output} ; done

fi; done

#................................................
#  Summary 19H1
#................................................

date ## Collect run summary for step 19H1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step19H1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tInput\tCandidates\tBest\tBase-file"
summary_output="${out_path_step22H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

fasta=`grep "CMD: minimap2" ${pbs_e} | cut -d" " -f13`;
main=`echo "$(basename ${fasta})" | sed 's/.\{6\}$//'`;

col01="${step}";                                                                                                    # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f5`;                                                                # Sample
col03=`grep -c "^>" ${fasta}`;                                                                                      # Input 
col04=`cut -f2 ${folder}/${main}/codpot_${main}/${main}.lncRNAclasses.txt | awk 'FNR>1{count++} END{print count}'`; # Candidates
col05=`cut -f2 ${folder}/${main}/codpot_${main}/${main}.lncRNAclasses.txt | sort | uniq | awk 'END{print NR}'`;     # Best
col06=`echo "$(basename ${pbs_e})"`;                                                                                # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}" >> ${summary_output} ; done

fi; done

#................................................
#  Summary 20H1
#................................................

date ## Collect run summary for step 20H1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step20H1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\t1model-EzLncPred\t2models-EzLncPred\t3models-EzLncPred\tMin2models-EzLncPred\tBest-FEELnc\tLncRNAs\tIntergenic\tBase-file"
summary_output="${out_path_step22H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${pbs_file} | while read pbs_p ; do

feelnc=`grep "grep -P" ${pbs_p} | cut -d" " -f4`;
main=`echo "$(basename ${feelnc})" | sed 's/.\{18\}$//'`;

col01="${step}";                                                                       # Step
col02=`echo "$(basename ${pbs_p})" | cut -d"_" -f5`;                                   # Sample
col03=`cut -f1 HYDRA_step21H1_*/${main}/${main}.CNCI.noncoding HYDRA_step21H1_*/${main}/${main}.CPAT.noncoding HYDRA_step21H1_*/${main}/${main}.CPC2.noncoding | sort | uniq -i -c | grep " 1 " | awk 'END{print NR}'`;                                             # 1model-EzLncPred
col04=`cut -f1 HYDRA_step21H1_*/${main}/${main}.CNCI.noncoding HYDRA_step21H1_*/${main}/${main}.CPAT.noncoding HYDRA_step21H1_*/${main}/${main}.CPC2.noncoding | sort | uniq -i -c | grep " 2 " | awk 'END{print NR}'`;                                             # 2models-EzLncPred 
col05=`cut -f1 HYDRA_step21H1_*/${main}/${main}.CNCI.noncoding HYDRA_step21H1_*/${main}/${main}.CPAT.noncoding HYDRA_step21H1_*/${main}/${main}.CPC2.noncoding | sort | uniq -i -c | grep " 3 " | awk 'END{print NR}'`;                                             # 3models-EzLncPred
col06=`awk 'END{print NR}' ${folder}/${main}/${main}.ezlncpred.min2-models.noncoding`; # Min2models-EzLncPred
col07=`awk 'END{print NR}' ${folder}/${main}/${main}.feelncmain.best`;                 # Best-FEELnc
col08=`awk 'END{print NR}' ${folder}/${main}/${main}.ezlncpred-and-feelnc.noncoding`;  # LncRNAs
col09=`awk 'END{print NR}' ${folder}/${main}/${main}.*.intergenic`;                    # Intergenic
col10=`echo "$(basename ${pbs_p})"`;                                                   # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}" >> ${summary_output} ; done

cut -f6,7,9- raw-rnaSPAdes_CDhit_identity98_cov90_merged_RPKM_lowcutoffest | sort | uniq -c

fi; done


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
sed -i 's,${bbduk},'"${bbduk}"',g' "$logfile"
sed -i 's,${memless},'"${memless}"',g' "$logfile"
sed -i 's,${cpuless},'"${cpuless}"',g' "$logfile"
sed -i 's,${ribosomal_rna_ref},'"${ribosomal_rna_ref}"',g' "$logfile"
sed -i 's,${out_path_step22H1_Bash},'"${out_path_step22H1_Bash}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v