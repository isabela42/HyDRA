#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Aug 10, 2023
Last modified on Jun 23, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 10H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step10H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 1 -c 1 -w "24:00:00"

## Input:

-i <path/to/input/files>    TSV file containing on:
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
#-->10 [H] Summary of quality control & read error correction steps

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
        i) input="${OPTARG}";;       # Input files including path
        p) pbs_stem="${OPTARG}";;    # Stem for PBS file names
        e) email="${OPTARG}";;       # Email for PBS job
        m) mem="${OPTARG}";;         # Memory required for PBS job
        c) ncpus="${OPTARG}";;       # Number of CPUS required for PBS job
        w) walltime="${OPTARG}";;    # Clock walltime required for PBS job
        h) Help ; exit;;             # Print Help and exit
        v) echo "${version}"; exit;; # Print version and exit
        ?) echo script usage: bash ipda_HYDRA_step10H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step10H1-to-pbs_${thislogdate}.txt

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
out_path_step10H1_Bash="HYDRA_step10H1_summary_Bash_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step10H1_Bash}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step10H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step10H1_Bash}"
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
echo "## Executing bash ipda_HYDRA_step10H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step10H1_Bash}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Summary 01L3
#................................................

date ## Collect run summary for step 01L3 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step01L3" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\t#Reads\tAveragePhred\tMedianPhred\tMinPhred\tMaxPhred\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

pbs=`echo ${pbs_o} | cut -d"." -f1`
inpfolder=`tail -n1 ${pbs}.pbs | cut -d" " -f5`

sample=`echo "$(basename ${pbs_o})" | cut -d"_" -f5-7`;

ls ${inpfolder}/${sample}*zip | cut -d"." -f1 | while read file; do 
base_file=`echo "$(basename ${file})"`; 
unzip -p ${file}.zip ${base_file}/fastqc_data.txt > fastqc_data.txt ;
lines=`grep -n ">>END_MODULE" fastqc_data.txt | cut -d":" -f1 | head -n2 | tail -n1`; y=$((lines-14)) ; z=$((y-1));
grep -A${y} ">>Per base sequence quality" fastqc_data.txt | tail -n${z} > temp.txt

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;                                            # Sample
col03=`grep "Total Sequences" fastqc_data.txt | cut -f2`;                                       # #Reads
col04=`cut -f2 temp.txt | awk '{s+=$1} END {print (NR ? s/NR : "NaN")}'`;                       # AveragePhred
col05=`cut -f2 temp.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`;                    # MedianPhred
col06=`cut -f2 temp.txt | awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }'`;  # MinPhred
col07=`cut -f2 temp.txt | awk '{if (NR == 1 || $1 > max) max = $1} END {print max}'`;           # MaxPhred
col08=`echo "$(basename ${pbs_o})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}" >> ${summary_output};

header="Step\tSample\tBase\tValues\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}_boxplot.tsv"
echo -e "${header}" > ${summary_output}

col01="${step}";                                          # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;      # Sample
col03=`cut -f1 temp.txt`;                                 # Base
col04=`cut -f2- temp.txt | tr '\t' '\n'`;                 # Values
col05=`echo "$(basename ${pbs_o})"`;                      # Base-file

echo ${col03} | tr ' ' '\n' | while read bp; do
echo ${col04} | tr ' ' '\n' | while read c4; do
echo -e "${col01}\t${col02}\t${bp}\t${c4}\t${col05}" >> ${summary_output};
done
done

done ; rm -rf fastqc_data.txt  temp.txt

done; 

fi; done

#................................................
#  Summary 01L5
#................................................

date ## Collect run summary for step 01L5 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step01L5" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tn_threshold\tTotal\tLongest\tSmallest\tn_plus\tup_t_n\t0_1k\t1k\t1k_2k\t2k_5k\t5k_10k\t10k_20k\t20k_30k\t30k_40k\t40k_50k\t50kplus\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

col01="${step}";                                        # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;    # Sample
col03=`tail -n16 ${pbs_o} | head -n1 | cut -f2`;        # n_threshold
col04=`tail -n15 ${pbs_o} | head -n1 | cut -f2`;        # Total
col05=`tail -n14 ${pbs_o} | head -n1 | cut -f2`;        # Longest
col06=`tail -n13 ${pbs_o} | head -n1 | cut -f2`;        # Smallest
col07=`tail -n12 ${pbs_o} | head -n1 | cut -f2`;        # n_plus
col08=`tail -n11 ${pbs_o} | head -n1 | cut -f2`;        # up_t_n
col09=`tail -n10 ${pbs_o} | head -n1 | cut -f2`;        # 0_1k
col10=`tail -n09 ${pbs_o} | head -n1 | cut -f2`;        # 1k
col11=`tail -n08 ${pbs_o} | head -n1 | cut -f2`;        # 1k_2k
col12=`tail -n07 ${pbs_o} | head -n1 | cut -f2`;        # 2k_5k
col13=`tail -n06 ${pbs_o} | head -n1 | cut -f2`;        # 5k_10k
col14=`tail -n05 ${pbs_o} | head -n1 | cut -f2`;        # 10k_20k
col15=`tail -n04 ${pbs_o} | head -n1 | cut -f2`;        # 20k_30k
col16=`tail -n03 ${pbs_o} | head -n1 | cut -f2`;        # 30k_40k
col17=`tail -n02 ${pbs_o} | head -n1 | cut -f2`;        # 40k_50k
col18=`tail -n01 ${pbs_o} | cut -f2`;                   # 50kplus
col19=`echo "$(basename ${pbs_o})"`;                    # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 01S1
#................................................

date ## Collect run summary for step 01S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step01S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tLane_Read\t#Reads\tAveragePhred\tMedianPhred\tMinPhred\tMaxPhred\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

pbs=`echo ${pbs_o} | cut -d"." -f1`
inpfolder=`tail -n1 ${pbs}.pbs | cut -d" " -f5`

sample=`echo "$(basename ${pbs_o})" | cut -d"_" -f5-8`;

ls ${inpfolder}/${sample}*zip | cut -d"." -f1 | while read file; do 
base_file=`echo "$(basename ${file})"`; 
unzip -p ${file}.zip ${base_file}/fastqc_data.txt > fastqc_data.txt ;
lines=`grep -n ">>END_MODULE" fastqc_data.txt | cut -d":" -f1 | head -n2 | tail -n1`; y=$((lines-14)) ; z=$((y-1));
grep -A${y} ">>Per base sequence quality" fastqc_data.txt | tail -n${z} > temp.txt

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f5`;                                            # Sample
col03=`grep "Filename" fastqc_data.txt | cut -f2 | cut -d"_" -f4-5 | cut -d"." -f1`;            # Lane_Read
col04=`grep "Total Sequences" fastqc_data.txt | cut -f2`;                                       # #Reads
col05=`cut -f2 temp.txt | awk '{s+=$1} END {print (NR ? s/NR : "NaN")}'`;                       # AveragePhred
col06=`cut -f2 temp.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`;                    # MedianPhred
col07=`cut -f2 temp.txt | awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }'`;  # MinPhred
col08=`cut -f2 temp.txt | awk '{if (NR == 1 || $1 > max) max = $1} END {print max}'`;           # MaxPhred
col09=`echo "$(basename ${pbs_o})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}" >> ${summary_output};

header="Step\tSample\tBase\tValues\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}_boxplot.tsv"
echo -e "${header}" > ${summary_output}

col01="${step}";                                          # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;      # Sample
col03=`cut -f1 temp.txt`;                                 # Base
col04=`cut -f2- temp.txt | tr '\t' '\n'`;                 # Values
col05=`echo "$(basename ${pbs_o})"`;                      # Base-file

echo ${col03} | tr ' ' '\n' | while read bp; do
echo ${col04} | tr ' ' '\n' | while read c4; do
echo -e "${col01}\t${col02}\t${bp}\t${c4}\t${col05}" >> ${summary_output};
done
done

done ; rm -rf fastqc_data.txt  temp.txt

done; 

fi; done

#................................................
#  Summary 02L1
#................................................

date ## Collect run summary for step 02L1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step02L1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tReads-trimmed-from-start\tPerc-reads-trimmed-from-start\tbp-removed-from-start\tReads-trimmed-from-end\tPerc-reads-trimmed-from-end\tbp-removed-from-end\tReads-split-based-middle-adapters\tIn-size\tOut-size\tIn-seq\tOut-seq\tIn-file\tOut-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

inpfile=`grep -A5 "## Run Porechop at" ${pbs_o} | tail -n1`
outfile=`tail -n2 ${pbs_o} | head -n1 | cut -d" " -f4`

col01="${step}"                                                                         # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;                                    # Sample
col03=`tail -n15 ${pbs_o} | head -n1 | cut -d" " -f3 | tr -d ','`;                      # In-reads
col04=`tail -n15 ${pbs_o} | head -n1 | cut -d" " -f1 | tr -d ','`;                      # Reads-trimmed-from-start
col05=`awk -v four="$col04" -v three="$col03" 'BEGIN {print ((four/three)*100)}'`;      # Perc-reads-trimmed-from-start
col06=`tail -n15 ${pbs_o} | head -n1 | cut -d" " -f11 | cut -d"(" -f2 | tr -d ','`;     # bp-removed-from-start
col07=`tail -n14 ${pbs_o} | head -n1 | cut -d" " -f1 | tr -d ','`;                      # Reads-trimmed-from-end
col08=`awk -v seven="$col07" -v three="$col03" 'BEGIN {print ((seven/three)*100)}'`;    # Perc-reads-trimmed-from-end
col09=`tail -n14 ${pbs_o} | head -n1 | cut -d" " -f11 | cut -d"(" -f2 | tr -d ','`;     # bp-removed-from-end
col10=`tail -n8 ${pbs_o} | head -n1 | cut -d" " -f1 | tr -d ','`;                       # Reads-split-based-middle-adapters
col11=`ls -lh ${inpfile} | cut -d" " -f5`;                                              # In-size
col12=`ls -lh ${outfile} | cut -d" " -f5`;                                              # Out-size
col13=`zcat -c ${inpfile} | awk 'END{print NR/4}'`;                                     # In-seq
col14=`zcat -c ${outfile} | awk 'END{print NR/4}'`;                                     # Out-seq
col15="${inpfile}";                                                                     # In-file
col16="${outfile}";                                                                     # Out-file
col17=`echo "$(basename ${pbs_o})"` ;                                                   # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 02L2
#................................................

date ## Collect run summary for step 02L2 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step02L2" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tOut-reads\tIn-size\tOut-size\tIn-seq\tOut-seq\tIn-file\tOut-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

pbs=`echo ${pbs_e} | cut -d"." -f1`
inpfile=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f3`
outfile=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f15`

col01="${step}"                                         # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f8`;    # Sample
col03=`tail -n1 ${pbs_e} | head -n1 | cut -d" " -f6`;   # In-reads
col04=`tail -n1 ${pbs_e} | head -n1 | cut -d" " -f2`;   # Out-reads
col05=`ls -lh ${inpfile} | cut -d" " -f5`;              # In-size
col06=`ls -lh ${outfile} | cut -d" " -f5`;              # Out-size
col07=`zcat -c ${inpfile} | awk 'END{print NR/4}'`;     # In-seq
col08=`zcat -c ${outfile} | awk 'END{print NR/4}'`;     # Out-seq
col09="${inpfile}";                                     # In-file
col10="${outfile}";                                     # Out-file
col11=`echo "$(basename ${pbs_e})"` ;                   # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 02L3
#................................................

date ## Collect run summary for step 02L3 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step02L3" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tIn-bases\tOut-reads\tPerc-out-reads\tOut-bases\tPerc-out-bases\tpolyA-trimmed-bp\tPerc-polyA-trimmed-bp\tIn-size\tOut-size\tIn-seq\tOut-seq\tIn-file\tOut-file\tBase-file"

summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

inpfile=`tail -n12 ${pbs_o} | head -n1 | cut -d" " -f5`
outfile=`tail -n12 ${pbs_o} | head -n1 | cut -d" " -f7`

col01="${step}"                                                                                     # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f9`;                                                # Sample
col03=`tail -n6 ${pbs_o} | head -n1 | cut -d" " -f17 | tr -d ','`;                                  # In-reads
col04=`tail -n3 ${pbs_o} | head -n1 | cut -d" " -f4 | tr -d ','`;                                   # In-bases
col05=`tail -n5 ${pbs_o} | head -n1 | cut -d" " -f8 | tr -d ','`;                                   # Out-reads
col06=`tail -n5 ${pbs_o} | head -n1 | cut -d" " -f9 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;   # Perc-out-reads
col07=`tail -n1 ${pbs_o} | cut -d" " -f5 | tr -d ','`;                                              # Out-bases
col08=`tail -n1 ${pbs_o} | cut -d" " -f7 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;              # Perc-out-bases
col09=`tail -n2 ${pbs_o} | head -n1 | cut -d" " -f15 | tr -d ','`;                                  # polyA-trimmed-bp
col10=`tail -n2 ${pbs_o} | head -n1 | cut -d" " -f17 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;  # Perc-polyA-trimmed-bp
col11=`ls -lh ${inpfile} | cut -d" " -f5`;                                                          # In-size
col12=`ls -lh ${outfile} | cut -d" " -f5`;                                                          # Out-size
col13=`zcat -c ${inpfile} | awk 'END{print NR/4}'`;                                                 # In-seq
col14=`zcat -c ${outfile} | awk 'END{print NR/4}'`;                                                 # Out-seq
col15="${inpfile}";                                                                                 # In-file
col16="${outfile}";                                                                                 # Out-file
col17=`echo "$(basename ${pbs_o})"` ;                                                               # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 02S1
#................................................

date ## Collect run summary for step 02S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step02S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tIn-bases\tCor-reads\tCor-bases\tPerc-cor-bases\tCor-bases-per-read\tKmers-stores\tWeak-kmer-threshold-rate\tBad-quality-threshold\tR1raw-size\tR1cor-size\tR2raw-size\tR2cor-size\tR1raw-seq\tR1cor-seq\tR2raw-seq\tR2cor-seq\tR1raw-file\tR1cor-file\tR2raw-file\tR2cor-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

R1raw=`grep -A1 "Put the kmers into bloom filter" ${pbs_e} | tail -n1 | cut -d" " -f14 | cut -d")" -f1`
R2raw=`grep -A1 "Put the kmers into bloom filter" ${pbs_e} | tail -n1 | cut -d" " -f17 | cut -d")" -f1`
R1cor=`grep -A1 "Put the kmers into bloom filter" ${pbs_e} | tail -n1 | cut -d" " -f14 | cut -d")" -f1 | cut -d"." -f1`
R2cor=`grep -A1 "Put the kmers into bloom filter" ${pbs_e} | tail -n1 | cut -d" " -f17 | cut -d")" -f1 | cut -d"." -f1`

col01="${step}";                                                                            # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f5` ;                                       # Sample
col03=`grep -B1 -Pe "^\tCorrected " ${pbs_e} | head -n1 | cut -d" " -f2`;                   # In-read
col04=`grep "Input:" ${pbs_e} | cut -f3 | cut -d" " -f1`;                                   # In-bases
col05=`grep "Input:" ${pbs_e} | cut -f2 | cut -d" " -f1`;                                   # Cor-reads
col06=`grep -Pe "^\tCorrected " ${pbs_e} | cut -d" " -f2`;                                  # Cor-bases
col07=`awk -v six="$col06" -v four="$col04" 'BEGIN {print ((six/four)*100)}'`;              # Perc-cor-bases
col08=`awk -v six="$col06" -v three="$col03" 'BEGIN {print ((six/three)*100)}'`;            # Cor-bases-per-read
col09=`grep -B4 -Pe "^\tCorrected " ${pbs_e} | head -n1 | cut -d" " -f2`;                   # Kmers-stores
col10=`grep -B3 -Pe "^\tCorrected " ${pbs_e} | head -n1 | cut -d" " -f5`;                   # Weak-kmer-threshold-rate
col11=`grep -B2 -Pe "^\tCorrected " ${pbs_e} | head -n1 | cut -d" " -f5 | cut -d"'" -f2`;   # Bad-quality-threshold
col12=`ls -lh ${R1raw} | cut -d" " -f5`;                                                    # R1raw-size (.fq.gz)
col13=`ls -lh ${R1cor}.cor.fq.gz | cut -d" " -f5`;                                          # R1cor-size (.fq.gz)
col14=`ls -lh ${R2raw} | cut -d" " -f5`;                                                    # R2raw-size (.fq.gz)
col15=`ls -lh ${R2cor}.cor.fq.gz | cut -d" " -f5`;                                          # R2cor-size (.fq.gz)
col16=`zcat -c ${R1raw} | awk 'END{print NR/4}'`;                                           # R1raw-seq
col17=`zcat -c ${R1cor}.cor.fq.gz | awk 'END{print NR/4}'`;                                 # R1cor-seq
col18=`zcat -c ${R2raw} | awk 'END{print NR/4}'`;                                           # R2raw-seq
col19=`zcat -c ${R2cor}.cor.fq.gz | awk 'END{print NR/4}'`;                                 # R2cor-seq
col20="${R1raw}";                                                                           # R1raw-file
col21="${R1cor}.cor.fq.gz";                                                                 # R1cor-file
col22="${R2raw}";                                                                           # R2raw-file
col23="${R2cor}.cor.fq.gz";                                                                 # R2cor-file
col24=`echo "$(basename ${pbs_e})"`;                                                        # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 02S2
#................................................

date ## Collect run summary for step 02S2 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step02S2" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-PE-reads\tRemoved-PE-reads\tOut-PE-reads\tR1-corrected\tR2-corrected\tPairs-corrected\tR1-unfixable\tR2-unfixable\tPairs-unfixable\tR1in-size\tR1out-size\tR2in-size\tR2out-size\tR1in-seq\tR1out-seq\tR2in-seq\tR2out-seq\tR1in-file\tR1out-file\tR2in-file\tR2out-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

col01="${step}";                                                                                                                        # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f5`;                                                                                    # Sample
col03=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n1 ${folder}/${file} | cut -d":" -f2 ; done`;             # In-PE-reads
col04=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n2 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # Removed-PE-reads
col05=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n3 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # Out-PE-reads
col06=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n4 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # R1-corrected
col07=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n5 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # R2-corrected
col08=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n6 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # Pairs-corrected
col09=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n7 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # R1-unfixable
col10=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n8 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # R2-unfixable
col11=`ls ${folder} | grep  "${col02}" | head -n1 | while read file; do head -n9 ${folder}/${file} | tail -n1 | cut -d":" -f2 ; done`;  # Pairs-unfixable

pbs=`echo ${pbs_e} | cut -d"." -f1`
R1inp=`tail -n6 ${pbs}.pbs | head -n1 | cut -d" " -f4`
R2inp=`tail -n6 ${pbs}.pbs | head -n1 | cut -d" " -f6`
R1out=`ls ${folder} | grep  "${col02}" | head -n2 | tail -n1`
R2out=`ls ${folder} | grep  "${col02}" | head -n3 | tail -n1`

col12=`ls -lh ${R1inp} | cut -d" " -f5`;                                                                                                # R1in-size
col13=`ls -lh ${folder}/${R1out} | cut -d" " -f5`;                                                                                      # R1out-size
col14=`ls -lh ${R2inp} | cut -d" " -f5`;                                                                                                # R2in-size
col15=`ls -lh ${folder}/${R2out} | cut -d" " -f5`;                                                                                      # R2out-size
col16=`zcat -c ${R1inp} | awk 'END{print NR/4}'`;                                                                                       # R1in-seq
col17=`awk 'END{print NR/4}' ${folder}/${R1out}`;                                                                                       # R1out-seq
col18=`zcat -c ${R2inp} | awk 'END{print NR/4}'`;                                                                                       # R2in-seq
col19=`awk 'END{print NR/4}' ${folder}/${R2out}`;                                                                                       # R2out-seq
col20="${R1inp}";                                                                                                                       # R1in-file
col21="${folder}/${R1out}";                                                                                                             # R1out-file
col22="${R2inp}";                                                                                                                       # R2in-file
col23="${folder}/${R2out}";                                                                                                             # R2out-file
col24=`echo "$(basename ${pbs_e})"`;                                                                                                    # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 03S1
#................................................

date ## Collect run summary for step 03S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step03S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tLane_Read\t#Reads\tAveragePhred\tMedianPhred\tMinPhred\tMaxPhred\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

pbs=`echo ${pbs_o} | cut -d"." -f1`
inpfolder=`tail -n1 ${pbs}.pbs | cut -d" " -f5`

sample=`echo "$(basename ${pbs_o})" | cut -d"_" -f6`;

ls ${inpfolder}/*${sample}*zip | cut -d"." -f1-2| while read file; do 
base_file=`echo "$(basename ${file})"`;
unzip -p ${file}.zip ${base_file}/fastqc_data.txt > fastqc_data.txt ;
lines=`grep -n ">>END_MODULE" fastqc_data.txt | cut -d":" -f1 | head -n2 | tail -n1`; y=$((lines-14)) ; z=$((y-1));
grep -A${y} ">>Per base sequence quality" fastqc_data.txt | tail -n${z} > temp.txt

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f6`;                                            # Sample
col03=`grep "Filename" fastqc_data.txt | cut -f2 | cut -d"_" -f4-5 | cut -d"." -f1`;            # Lane_Read
col04=`grep "Total Sequences" fastqc_data.txt | cut -f2`;                                       # #Reads
col05=`cut -f2 temp.txt | awk '{s+=$1} END {print (NR ? s/NR : "NaN")}'`;                       # AveragePhred
col06=`cut -f2 temp.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`;                    # MedianPhred
col07=`cut -f2 temp.txt | awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }'`;  # MinPhred
col08=`cut -f2 temp.txt | awk '{if (NR == 1 || $1 > max) max = $1} END {print max}'`;           # MaxPhred
col09=`echo "$(basename ${pbs_o})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}" >> ${summary_output};

done ; rm -rf fastqc_data.txt temp.txt

done; 

fi; done

#................................................
#  Summary 03L3
#................................................

date ## Collect run summary for step 03L3 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step03L3" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\t#Reads\tAveragePhred\tMedianPhred\tMinPhred\tMaxPhred\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

pbs=`echo ${pbs_o} | cut -d"." -f1`
inpfolder=`tail -n1 ${pbs}.pbs | cut -d" " -f5`

sample=`echo "$(basename ${pbs_o})" | cut -d"_" -f8-10`;

ls ${inpfolder}/*${sample}*zip | cut -d"." -f1 | while read file; do 
base_file=`echo "$(basename ${file})"`; 
unzip -p ${file}.zip ${base_file}/fastqc_data.txt > fastqc_data.txt ;
lines=`grep -n ">>END_MODULE" fastqc_data.txt | cut -d":" -f1 | head -n2 | tail -n1`; y=$((lines-14)) ; z=$((y-1));
grep -A${y} ">>Per base sequence quality" fastqc_data.txt | tail -n${z} > temp.txt

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f10`;                                           # Sample
col03=`grep "Total Sequences" fastqc_data.txt | cut -f2`;                                       # #Reads
col04=`cut -f2 temp.txt | awk '{s+=$1} END {print (NR ? s/NR : "NaN")}'`;                       # AveragePhred
col05=`cut -f2 temp.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`;                    # MedianPhred
col06=`cut -f2 temp.txt | awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }'`;  # MinPhred
col07=`cut -f2 temp.txt | awk '{if (NR == 1 || $1 > max) max = $1} END {print max}'`;           # MaxPhred
col08=`echo "$(basename ${pbs_o})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}" >> ${summary_output};

done ; rm -rf fastqc_data.txt  temp.txt

done; 

fi; done

#................................................
#  Summary 03L5
#................................................

date ## Collect run summary for step 03L5 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step03L5" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tn_threshold\tTotal\tLongest\tSmallest\tn_plus\tup_t_n\t0_1k\t1k\t1k_2k\t2k_5k\t5k_10k\t10k_20k\t20k_30k\t30k_40k\t40k_50k\t50kplus\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

col01="${step}";                                        # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f10`;   # Sample
col03=`tail -n16 ${pbs_o} | head -n1 | cut -f2`;        # n_threshold
col04=`tail -n15 ${pbs_o} | head -n1 | cut -f2`;        # Total
col05=`tail -n14 ${pbs_o} | head -n1 | cut -f2`;        # Longest
col06=`tail -n13 ${pbs_o} | head -n1 | cut -f2`;        # Smallest
col07=`tail -n12 ${pbs_o} | head -n1 | cut -f2`;        # n_plus
col08=`tail -n11 ${pbs_o} | head -n1 | cut -f2`;        # up_t_n
col09=`tail -n10 ${pbs_o} | head -n1 | cut -f2`;        # 0_1k
col10=`tail -n09 ${pbs_o} | head -n1 | cut -f2`;        # 1k
col11=`tail -n08 ${pbs_o} | head -n1 | cut -f2`;        # 1k_2k
col12=`tail -n07 ${pbs_o} | head -n1 | cut -f2`;        # 2k_5k
col13=`tail -n06 ${pbs_o} | head -n1 | cut -f2`;        # 5k_10k
col14=`tail -n05 ${pbs_o} | head -n1 | cut -f2`;        # 10k_20k
col15=`tail -n04 ${pbs_o} | head -n1 | cut -f2`;        # 20k_30k
col16=`tail -n03 ${pbs_o} | head -n1 | cut -f2`;        # 30k_40k
col17=`tail -n02 ${pbs_o} | head -n1 | cut -f2`;        # 40k_50k
col18=`tail -n01 ${pbs_o} | cut -f2`;                   # 50kplus
col19=`echo "$(basename ${pbs_o})"`;                    # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 04L1
#................................................

date ## Collect run summary for step 04L1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step04L1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\t$\tA\tC\tG\tT\tN\t1R1-size\t1R2-size\t2R1-size\t2R2-size\t3R1-size\t3R2-size\t4R1-size\t4R2-size\t5R1-size\t5R2-size\tOut-size\t1R1-file\t1R2-file\t2R1-file\t2R2-file\t3R1-file\t3R2-file\t4R1-file\t4R2-file\t5R1-file\t5R2-file\tOut-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

pbs=`echo ${pbs_e} | cut -d"." -f1`
oneR1file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f3`
oneR2file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f4`
twoR1file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f5`
twoR2file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f6`
threeR1file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f7`
threeR2file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f8`
fourR1file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f9`
fourR2file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f10`
fiveR1file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f11`
fiveR2file=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f12`
outfile=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f33`

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f4`;                                            # Sample
col03=`tail -n4 ${pbs_e} | head -n1 | cut -d"(" -f3 | cut -d"," -f1`;                           # $
col04=`tail -n4 ${pbs_e} | head -n1 | cut -d"(" -f3 | cut -d"," -f2 | tr -d ' '`;               # A
col05=`tail -n4 ${pbs_e} | head -n1 | cut -d"(" -f3 | cut -d"," -f3 | tr -d ' '`;               # C
col06=`tail -n4 ${pbs_e} | head -n1 | cut -d"(" -f3 | cut -d"," -f4 | tr -d ' '`;               # G
col07=`tail -n4 ${pbs_e} | head -n1 | cut -d"(" -f3 | cut -d"," -f5 | tr -d ' '`;               # T
col08=`tail -n4 ${pbs_e} | head -n1 | cut -d"(" -f3 | cut -d"," -f6 | tr -d ' ' | tr -d ')'`;   # N
col09=`ls -lh ${oneR1file} | cut -d" " -f5`;                                                    # 1R1-size
col10=`ls -lh ${oneR2file} | cut -d" " -f5`;                                                    # 1R2-size
col11=`ls -lh ${twoR1file} | cut -d" " -f5`;                                                    # 2R1-size
col12=`ls -lh ${twoR2file} | cut -d" " -f5`;                                                    # 2R2-size
col13=`ls -lh ${threeR1file} | cut -d" " -f5`;                                                  # 3R1-size
col14=`ls -lh ${threeR2file} | cut -d" " -f5`;                                                  # 3R2-size
col15=`ls -lh ${fourR1file} | cut -d" " -f5`;                                                   # 4R1-size
col16=`ls -lh ${fourR2file} | cut -d" " -f5`;                                                   # 4R2-size
col17=`ls -lh ${fiveR1file} | cut -d" " -f5`;                                                   # 5R1-size
col18=`ls -lh ${fiveR2file} | cut -d" " -f5`;                                                   # 5R2-size
col19=`ls -lh ${outfile} | cut -d" " -f5`;                                                      # Out-size
col20="${oneR1file}";                                                                           # 1R1-file
col21="${oneR2file}";                                                                           # 1R2-file
col22="${twoR1file}";                                                                           # 2R1-file
col23="${twoR2file}";                                                                           # 2R2-file
col24="${threeR1file}";                                                                         # 3R1-file
col25="${threeR2file}";                                                                         # 3R2-file
col26="${fourR1file}";                                                                          # 4R1-file
col27="${fourR2file}";                                                                          # 4R2-file
col28="${fiveR1file}";                                                                          # 5R1-file
col29="${fiveR2file}";                                                                          # 5R2-file
col30="${outfile}";                                                                             # Out-file
col31=`echo "$(basename ${pbs_e})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}\t${col25}\t${col26}\t${col27}\t${col28}\t${col29}\t${col30}\t${col31}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 04L2
#................................................

date ## Collect run summary for step 04L2 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step04L2" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\t$\tA\tC\tG\tT\tN\tRLE-BWT-byte-length\tStatus\tIn-size\tOut-size\tIn-file\tOut-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

pbs=`echo ${pbs_e} | cut -d"." -f1`
inpfile=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f3`
outfile=`tail -n4 ${pbs}.pbs | head -n1 | cut -d" " -f4`

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f5`;                                            # Sample
col03=`tail -n3 ${pbs_e} | head -n1 | cut -d"[" -f3 | cut -d"," -f1`;                           # $
col04=`tail -n3 ${pbs_e} | head -n1 | cut -d"[" -f3 | cut -d"," -f2 | tr -d ' '`;               # A
col05=`tail -n3 ${pbs_e} | head -n1 | cut -d"[" -f3 | cut -d"," -f3 | tr -d ' '`;               # C
col06=`tail -n3 ${pbs_e} | head -n1 | cut -d"[" -f3 | cut -d"," -f4 | tr -d ' '`;               # G
col07=`tail -n3 ${pbs_e} | head -n1 | cut -d"[" -f3 | cut -d"," -f5 | tr -d ' '`;               # T
col08=`tail -n3 ${pbs_e} | head -n1 | cut -d"[" -f3 | cut -d"," -f6 | tr -d ' ' | tr -d ']'`;   # N
col09=`tail -n2 ${pbs_e} | head -n1 | cut -d":" -f6 | tr -d ' '`;                               # RLE-BWT-byte-length
col10=`tail -n1 ${pbs_e} | cut -d" " -f5-`;                                                     # Status
col11=`ls -lh ${inpfile} | cut -d" " -f5`;                                                      # In-size
col12=`ls -lh ${outfile} | cut -d" " -f5`;                                                      # Out-size
col13="${inpfile}";                                                                             # In-file
col14="${outfile}";                                                                             # Out-file
col15=`echo "$(basename ${pbs_e})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 04L3
#................................................

date ## Collect run summary for step 04L3 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step04L3" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tTotal-reads\tRange\tIn-npy-size\tIn-reads-size\tOut-size\tIn-reads-seq\tOut-seq\tIn-npy-file\tIn-reads-file\tOut-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

inpnpy=`head -n2 ${pbs_e}| tail -n1 | cut -d"\"" -f2`
inpreads=`head -n3 ${pbs_e}| tail -n1 | cut -d"\"" -f2`
outfile=`head -n4 ${pbs_e}| tail -n1 | cut -d"\"" -f2`

col01="${step}";                                       # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f10`;  # Sample
col03=`tail -n1 ${pbs_e} | cut -d" " -f7`;             # Total-reads
col04=`tail -n1 ${pbs_e} | cut -d" " -f12-`;           # Range
col05=`ls -lh ${inpnpy} | cut -d" " -f5`;              # In-npy-size
col06=`ls -lh ${inpreads} | cut -d" " -f5`;            # In-reads-size
col07=`ls -lh ${outfile}.gz | cut -d" " -f5`;          # Out-size
col08=`zcat -c ${inpreads} | awk 'END{print NR/4}'`;   # In-reads-seq
col09=`zcat -c ${outfile}.gz | awk 'END{print NR/2}'`; # Out-seq
col10="${inpnpy}";                                     # In-npy-file
col11="${inpreads}";                                   # In-reads-file
col12="${outfile}.gz";                                 # Out-file
col13=`echo "$(basename ${pbs_e})"`;                   # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 04S1
#................................................

date ## Collect run summary for step 04S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step04S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tPrefix-pairs\tForward-reverse-seqs\tForward-only-seqs\tReverse-only-seqs\tInput-read-pairs\tBoth-surviving\tPerc-both-surviving\tForward-only-surviving\tPerc-forward-only-surviving\tReverse-only-surviving\tPerc-reverse-only-surviving\tDropped\tPerc-dropped\tStatus\tR1in-size\tR1out-size\tR2in-size\tR2out-size\tR1in-seq\tR1out-seq\tR2in-seq\tR2out-seq\tR1in-file\tR1out-file\tR2in-file\tR2out-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

R1inp=`grep -A1 "TrimmomaticPE: Started with arguments:" ${pbs_e} | tail -n1 | cut -d" " -f5`
R2inp=`grep -A1 "TrimmomaticPE: Started with arguments:" ${pbs_e} | tail -n1 | cut -d" " -f6`
R1out=`grep -A1 "TrimmomaticPE: Started with arguments:" ${pbs_e} | tail -n1 | cut -d" " -f7`
R2out=`grep -A1 "TrimmomaticPE: Started with arguments:" ${pbs_e} | tail -n1 | cut -d" " -f9`

col01="${step}"                                                                                                                 # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f6`;                                                                            # Sample
col03=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n1 | cut -d" " -f3`;                                                          # Prefix-pairs
col04=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n1 | cut -d" " -f6`;                                                          # Forward-reverse-seqs
col05=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n1 | cut -d" " -f9`;                                                          # Forward-only-seqs
col06=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n1 | cut -d" " -f13`;                                                         # Reverse-only-seqs
col07=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f4`;                                               # Input-read-pairs
col08=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f7`;                                               # Both-surviving
col09=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f8 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;   # Perc-both-surviving
col10=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f12`;                                              # Forward-only-surviving
col11=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f13 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;  # Perc-forward-only-surviving
col12=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f17`;                                              # Reverse-only-surviving
col13=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f18 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;  # Perc-reverse-only-surviving
col14=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f20`;                                              # Dropped
col15=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n2 | tail -n1 | cut -d" " -f21 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;  # Perc-dropped
col16=`grep -A2 "^ILLUMINACLIP:" ${pbs_e} | head -n3 | tail -n1`;                                                               # Status
col17=`ls -lh ${R1inp} | cut -d" " -f5`;                                                                                        # R1in-size
col18=`ls -lh ${R1out} | cut -d" " -f5`;                                                                                        # R1out-size
col19=`ls -lh ${R2inp} | cut -d" " -f5`;                                                                                        # R2in-size
col20=`ls -lh ${R2out} | cut -d" " -f5`;                                                                                        # R2out-size
col21=`awk 'END{print NR/4}' ${R1inp}`;                                                                                         # R1in-seq
col22=`zcat -c ${R1out} | awk 'END{print NR/4}'`;                                                                               # R1out-seq
col23=`awk 'END{print NR/4}' ${R2inp}`;                                                                                         # R2in-seq
col24=`zcat -c ${R2out} | awk 'END{print NR/4}'`;                                                                               # R2out-seq
col25="${R1inp}";                                                                                                               # R1in-file
col26="${R1out}";                                                                                                               # R1out-file
col27="${R2inp}";                                                                                                               # R2in-file
col28="${R2out}";                                                                                                               # R2out-file
col29=`echo "$(basename ${pbs_e})"`;                                                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}\t${col25}\t${col26}\t${col27}\t${col28}\t${col29}" >> ${summary_output}; done

echo -e "Count\tSum\tAverage\tMedian\tSmaller\tBigger" >> ${out_path_step10H1_Bash}/summary_${run_name}_statistics.tsv
cut -f9 ${summary_output} | grep -v "Perc-both-surviving" | tr -d '%' | sort -n | awk 'BEGIN { c = 0; sum = 0; }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; }
  END { ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ]; } else { median = ( a[c/2] + a[c/2-1] ) / 2; }
  OFS="\t"; print c"\t"sum"\t"ave"\t"median"\t"a[0]"\t"a[c-1]; }' >> ${out_path_step10H1_Bash}/summary_${run_name}_statistics.tsv

fi; done

#................................................
#  Summary 04S2
#................................................

date ## Collect run summary for step 04S2 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step04S2" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tIn-bases\tQ-trimmed-reads\tPerc-Q-trimmed-reads\tQ-Trimmed-bases\tPerc-Q-Trimmed-bases\tLow-Q-discards-reads\tPerc-low-Q-discards-reads\tLow-Q-discards-bases\tPerc-low-Q-discards-bases\tTotal-removed-reads\tPerc-total-removed-reads\tTotal-removed-bases\tPerc-total-removed-bases\tOut-reads\tPerc-out-reads\tOut-bases\tPerc-out-bases\tR1in-size\tR1failed-size\tR2in-size\tR2failed-size\tR1in-seq\tR1failed-seq\tR2in-seq\tR2failed-seq\tR1in-file\tR1failed-file\tR2in-file\tR2failed-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

R1inp=`grep -B13 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f5 | cut -b4- | cut -d"," -f1`
R2inp=`grep -B13 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f6 | cut -b5- | cut -d"," -f1`
R1out=`grep -B13 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f12 | cut -b6- | cut -d"," -f1`
R2out=`grep -B13 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f13 | cut -b7- | cut -d"," -f1`

col01="${step}";                                                                                                                                                # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f7`;                                                                                                            # Sample
col03=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n5 | tail -n1 | cut -f2 | cut -d" " -f1`;                                                # In-reads
col04=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n5 | tail -n1 | cut -f4 | cut -d" " -f1`;                                                # In-bases
col05=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f2 | cut -d" " -f1`;                                                # Q-trimmed-reads
col06=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-Q-trimmed-reads
col07=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f3 | cut -d" " -f1`;                                                # Q-Trimmed-bases
col08=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-Q-Trimmed-bases
col09=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f2 | cut -d" " -f1`;                                                # Low-Q-discards-reads
col10=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1`;                # Perc-low-Q-discards-reads
col11=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f3 | cut -d" " -f1`;                                                # Low-Q-discards-bases
col12=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1`;                # Perc-low-Q-discards-bases
col13=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f2 | cut -d" " -f1`;                                                # Total-removed-reads
col14=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-total-removed-reads
col15=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f3 | cut -d" " -f1`;                                                # Total-removed-bases
col16=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-total-removed-bases
col17=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n9 | tail -n1 | cut -f2 | cut -d" " -f1`;                                                # Out-reads
col18=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n9 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-out-reads
col19=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n9 | tail -n1 | cut -f3 | cut -d" " -f1`;                                                # Out-bases
col20=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n9 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-out-bases
col21=`ls -lh ${R1inp} | cut -d" " -f5`;                                                                                                                        # R1in-size
col22=`ls -lh ${R1out} | cut -d" " -f5`;                                                                                                                        # R1failed-size
col23=`ls -lh ${R2inp} | cut -d" " -f5`;                                                                                                                        # R2in-size
col24=`ls -lh ${R2out} | cut -d" " -f5`;                                                                                                                        # R2failed-size
col25=`zcat -c ${R1inp} | awk 'END{print NR/4}'`;                                                                                                               # R1in-seq
col26=`zcat -c ${R1out} | awk 'END{print NR/4}'`;                                                                                                               # R1failed-seq
col27=`zcat -c ${R2inp} | awk 'END{print NR/4}'`;                                                                                                               # R2in-seq
col28=`zcat -c ${R2out} | awk 'END{print NR/4}'`;                                                                                                               # R2failed-seq
col29="${R1inp}";                                                                                                                                               # R1in-file
col30="${R1out}";                                                                                                                                               # R1failed-file
col31="${R2inp}";                                                                                                                                               # R2in-file
col32="${R2out}";                                                                                                                                               # R2failed-file
col33=`echo "$(basename ${pbs_e})"`;                                                                                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}\t${col25}\t${col26}\t${col27}\t${col28}\t${col29}\t${col30}\t${col31}\t${col32}\t${col33}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 05L1
#................................................

date ## Collect run summary for step 05L1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step05L1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tn_threshold\tTotal\tLongest\tSmallest\tn_plus\tup_t_n\t0_1k\t1k\t1k_2k\t2k_5k\t5k_10k\t10k_20k\t20k_30k\t30k_40k\t40k_50k\t50kplus\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

col01="${step}";                                        # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f11`;   # Sample
col03=`tail -n19 ${pbs_o} | head -n1 | cut -f2`;        # n_threshold
col04=`tail -n18 ${pbs_o} | head -n1 | cut -f2`;        # Total
col05=`tail -n17 ${pbs_o} | head -n1 | cut -f2`;        # Longest
col06=`tail -n16 ${pbs_o} | head -n1 | cut -f2`;        # Smallest
col07=`tail -n15 ${pbs_o} | head -n1 | cut -f2`;        # n_plus
col08=`tail -n14 ${pbs_o} | head -n1 | cut -f2`;        # up_t_n
col09=`tail -n13 ${pbs_o} | head -n1 | cut -f2`;        # 0_1k
col10=`tail -n12 ${pbs_o} | head -n1 | cut -f2`;        # 1k
col11=`tail -n11 ${pbs_o} | head -n1 | cut -f2`;        # 1k_2k
col12=`tail -n10 ${pbs_o} | head -n1 | cut -f2`;        # 2k_5k
col13=`tail -n09 ${pbs_o} | head -n1 | cut -f2`;        # 5k_10k
col14=`tail -n08 ${pbs_o} | head -n1 | cut -f2`;        # 10k_20k
col15=`tail -n07 ${pbs_o} | head -n1 | cut -f2`;        # 20k_30k
col16=`tail -n06 ${pbs_o} | head -n1 | cut -f2`;        # 30k_40k
col17=`tail -n05 ${pbs_o} | head -n1 | cut -f2`;        # 40k_50k
col18=`tail -n04 ${pbs_o} | head -n1 | cut -f2`;        # 50kplus
col19=`echo "$(basename ${pbs_o})"`;                    # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 05S1
#................................................

date ## Collect run summary for step 05S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step05S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tLane_Read\t#Reads\tAveragePhred\tMedianPhred\tMinPhred\tMaxPhred\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

pbs=`echo ${pbs_o} | cut -d"." -f1`
inpfolder=`tail -n1 ${pbs}.pbs | cut -d" " -f5`

sample=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;

ls ${inpfolder}/*${sample}*zip | cut -d"." -f1| while read file; do 
base_file=`echo "$(basename ${file})"`;
unzip -p ${file}.zip ${base_file}/fastqc_data.txt > fastqc_data.txt ;
lines=`grep -n ">>END_MODULE" fastqc_data.txt | cut -d":" -f1 | head -n2 | tail -n1`; y=$((lines-14)) ; z=$((y-1));
grep -A${y} ">>Per base sequence quality" fastqc_data.txt | tail -n${z} > temp.txt

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f7`;                                            # Sample
col03=`grep "Filename" fastqc_data.txt | cut -f2 | cut -d"_" -f6 | cut -d"." -f1`;              # Lane_Read
col04=`grep "Total Sequences" fastqc_data.txt | cut -f2`;                                       # #Reads
col05=`cut -f2 temp.txt | awk '{s+=$1} END {print (NR ? s/NR : "NaN")}'`;                       # AveragePhred
col06=`cut -f2 temp.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`;                    # MedianPhred
col07=`cut -f2 temp.txt | awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }'`;  # MinPhred
col08=`cut -f2 temp.txt | awk '{if (NR == 1 || $1 > max) max = $1} END {print max}'`;           # MaxPhred
col09=`echo "$(basename ${pbs_o})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}" >> ${summary_output};

done ; rm -rf fastqc_data.txt temp.txt

done; 

fi; done

#................................................
#  Summary 06L1
#................................................

date ## Collect run summary for step 06L1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step06L1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tIn-bases\trRNA-reads\tPerc-rRNA-reads\trRNA-bases\tPerc-rRNA-bases\tTotal-removed-reads\tPerc-total-removed-reads\tTotal-removed-bases\tPerc-total-removed-bases\tOut-reads\tPerc-out-reads\tOut-bases\tPerc-out-bases\tIn-size\trRNA-size\tRibodepleted-size\tIn-seq\trRNA-seq\tRibodepleted-seq\tIn-file\trRNA-file\tRibodepleted-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

inp=`grep -B11 "^Input is being processed as unpaired" ${pbs_e} | head -n1 | cut -d" " -f5 | cut -b4- | cut -d"," -f1`
rna=`grep -B11 "^Input is being processed as unpaired" ${pbs_e} | head -n1 | cut -d" " -f6 | cut -b6- | cut -d"," -f1`
dep=`grep -B11 "^Input is being processed as unpaired" ${pbs_e} | head -n1 | cut -d" " -f7 | cut -b5- | cut -d"," -f1`

col01="${step}";                                                                                                                                     # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f11`;                                                                                                 # Sample
col03=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n8 | tail -n1 | cut -f2 | cut -d" " -f1`;                                  # In-reads
col04=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n8 | tail -n1 | cut -f4 | cut -d" " -f1`;                                  # In-bases
col05=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n9 | tail -n1 | cut -f2 | cut -d" " -f1`;                                  # rRNA-reads
col06=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n9 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;  # Perc-rRNA-reads
col07=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n9 | tail -n1 | cut -f3 | cut -d" " -f1`;                                  # rRNA-bases
col08=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n9 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;  # Perc-rRNA-bases
col09=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n10 | tail -n1 | cut -f2 | cut -d" " -f1`;                                 # Total-removed-reads
col10=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n10 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`; # Perc-total-removed-read
col11=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n10 | tail -n1 | cut -f3 | cut -d" " -f1`;                                 # Total-removed-bases
col12=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n10 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`; # Perc-total-removed-bp
col13=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n11 | tail -n1 | cut -f2 | cut -d" " -f1`;                                 # Out-reads
col14=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n11 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`; # Perc-out-reads
col15=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n11 | tail -n1 | cut -f3 | cut -d" " -f1`;                                 # Out-bases
col16=`grep -A10 "^Input is being processed as unpaired" ${pbs_e} | head -n11 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`; # Perc-out-bases
col17=`ls -lh ${inp} | cut -d" " -f5`;                                                                                                               # In-size
col18=`ls -lh ${rna} | cut -d" " -f5`;                                                                                                               # rRNA-size
col19=`ls -lh ${dep} | cut -d" " -f5`;                                                                                                               # Ribodepleted-size
col20=`zcat -c ${inp} | awk 'END{print NR/2}'`;                                                                                                      # In-seq
col21=`awk 'END{print NR/2}' ${rna}`;                                                                                                                # rRNA-seq
col22=`awk 'END{print NR/2}' ${dep}`;                                                                                                                # Ribodepleted-seq
col23="${inp}";                                                                                                                                      # In-file
col24="${rna}";                                                                                                                                      # rRNA-file
col25="${dep}";                                                                                                                                      # Ribodepleted-size
col26=`echo "$(basename ${pbs_e})"`;                                                                                                                 # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}\t${col25}\t${col26}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 06S1
#................................................

date ## Collect run summary for step 06S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step06S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tIn-reads\tIn-bases\trRNA-reads\tPerc-rRNA-reads\trRNA-bases\tPerc-rRNA-bases\tTotal-removed-reads\tPerc-total-removed-reads\tTotal-removed-bases\tPerc-total-removed-bases\tOut-reads\tPerc-out-reads\tOut-bases\tPerc-out-bases\tR1in-size\tR1rRNA-size\tR1ribodepleted-size\tR2in-size\tR2rRNA-size\tR2ribodepleted-size\tR1in-seq\tR1rRNA-seq\tR1ribodepleted-seq\tR2in-seq\tR2rRNA-seq\tR2ribodepleted-seq\tR1in-file\tR1rRNA-file\tR1ribodepleted-file\tR2in-file\tR2rRNA-file\tR2ribodepleted-file\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

R1inp=`grep -B11 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f5 | cut -b4- | cut -d"," -f1`
R2inp=`grep -B11 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f6 | cut -b5- | cut -d"," -f1`
R1rna=`grep -B11 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f7 | cut -b6- | cut -d"," -f1`
R2rna=`grep -B11 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f8 | cut -b7- | cut -d"," -f1`
R1dep=`grep -B11 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f9 | cut -b5- | cut -d"," -f1`
R2dep=`grep -B11 "^Input is being processed as paired" ${pbs_e} | head -n1 | cut -d" " -f10 | cut -b6- | cut -d"," -f1`

col01="${step}";                                                                                                                                    # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f7`;                                                                                                # Sample
col03=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n5 | tail -n1 | cut -f2 | cut -d" " -f1`;                                    # In-reads
col04=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n5 | tail -n1 | cut -f4 | cut -d" " -f1`;                                    # In-bases
col05=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f2 | cut -d" " -f1`;                                    # rRNA-reads
col06=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-rRNA-reads
col07=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f3 | cut -d" " -f1`;                                    # rRNA-bases
col08=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n6 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-rRNA-bases
col09=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f2 | cut -d" " -f1`;                                    # Total-removed-reads
col10=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-total-removed-reads
col11=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f3 | cut -d" " -f1`;                                    # Total-removed-bases
col12=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n7 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-total-removed-bases
col13=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f2 | cut -d" " -f1`;                                    # Out-reads
col14=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f2 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-out-reads
col15=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f3 | cut -d" " -f1`;                                    # Out-bases
col16=`grep -A8 "^Input is being processed as paired" ${pbs_e} | head -n8 | tail -n1 | cut -f3 | cut -d" " -f3 | cut -d"(" -f2 | cut -d")" -f1 | tr -d '%'`;    # Perc-out-bases
col17=`ls -lh ${R1inp} | cut -d" " -f5`;                                                                                                            # R1in-size
col18=`ls -lh ${R1rna} | cut -d" " -f5`;                                                                                                            # R1rRNA-size
col19=`ls -lh ${R1dep} | cut -d" " -f5`;                                                                                                            # R1ribodepleted-size
col20=`ls -lh ${R2inp} | cut -d" " -f5`;                                                                                                            # R2in-size
col21=`ls -lh ${R2rna} | cut -d" " -f5`;                                                                                                            # R2rRNA-size
col22=`ls -lh ${R2dep} | cut -d" " -f5`;                                                                                                            # R2ribodepleted-size
col23=`zcat -c ${R1inp} | awk 'END{print NR/4}'`;                                                                                                   # R1in-seq
col24=`zcat -c ${R1rna} | awk 'END{print NR/4}'`;                                                                                                   # R1rRNA-seq
col25=`zcat -c ${R1dep} | awk 'END{print NR/4}'`;                                                                                                   # R1ribodepleted-seq
col26=`zcat -c ${R2inp} | awk 'END{print NR/4}'`;                                                                                                   # R2in-seq
col27=`zcat -c ${R2rna} | awk 'END{print NR/4}'`;                                                                                                   # R2rRNA-seq
col28=`zcat -c ${R2dep} | awk 'END{print NR/4}'`;                                                                                                   # R2ribodepleted-seq
col29="${R1inp}";                                                                                                                                   # R1in-file
col30="${R1rna}";                                                                                                                                   # R1rRNA-file
col31="${R1dep}";                                                                                                                                   # R1ribodepleted-size
col32="${R2inp}";                                                                                                                                   # R2in-file
col33="${R2rna}";                                                                                                                                   # R2rRNA-file
col34="${R2dep}";                                                                                                                                   # R2ribodepleted-size
col35=`echo "$(basename ${pbs_e})"`;                                                                                                                # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}\t${col24}\t${col25}\t${col26}\t${col27}\t${col28}\t${col29}\t${col30}\t${col31}\t${col32}\t${col33}\t${col34}\t${col35}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 07L1
#................................................

date ## Collect run summary for step 07L1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step07L1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tTotal-reads\tOveral-alignment-rate\tTotal-UP\tPerc-total-UP\tTotal-UP-aligned-0x\tPerc-total-UP-aligned-0x\tTotal-UP-aligned-1x\tPerc-total-UP-aligned-1x\tTotal-UP-aligned->1x\tPerc-total-UP-aligned->1x\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

col01="${step}";                                                                             # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f12`;                                        # Sample
col03=`tail -n7 ${pbs_e} | head -n1 | cut -d" " -f1`;                                        # Total-reads
col04=`tail -n2 ${pbs_e} | head -n1 | cut -d" " -f1 | tr -d '%'`;                            # Overal-alignment-rate
col05=`tail -n6 ${pbs_e} | head -n1 | cut -d" " -f3`;                                        # Total-UP
col06=`tail -n6 ${pbs_e} | head -n1 | cut -d" " -f4 | tr -d '(' | tr -d ')' | tr -d '%'`;    # Perc-total-UP
col07=`tail -n5 ${pbs_e} | head -n1 | cut -d" " -f5`;                                        # Total-UP-aligned-0x
col08=`tail -n5 ${pbs_e} | head -n1 | cut -d" " -f6 | tr -d '(' | tr -d ')' | tr -d '%'`;    # Perc-total-UP-aligned-0x
col09=`tail -n4 ${pbs_e} | head -n1 | cut -d" " -f5`;                                        # Total-UP-aligned-1x
col10=`tail -n4 ${pbs_e} | head -n1 | cut -d" " -f6 | tr -d '(' | tr -d ')' | tr -d '%'`;    # Perc-total-UP-aligned-1x
col11=`tail -n3 ${pbs_e} | head -n1 | cut -d" " -f5`;                                        # Total-UP-aligned->1x
col12=`tail -n3 ${pbs_e} | head -n1 | cut -d" " -f6 | tr -d '(' | tr -d ')' | tr -d '%'`;    # Perc-total-UP-aligned->1x
col13=`echo "$(basename ${pbs_e})"`;                                                         # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 07S1
#................................................

date ## Collect run summary for step 07S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step07S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tTotal-reads\tOveral-alignment-rate\tTotal-PE\tPerc-total-PE\tTotal-PE-aligned-concordantly-0x\tPerc-total-PE-aligned-concordantly-0x\tTotal-PE-aligned-concordantly-0x-aligned-discordantly-1x\tPerc-total-PE-aligned-concordantly-0x-aligned-discordantly-1x\tTotal-PE-aligned-concordantly-1x\tPerc-total-PE-aligned-concordantly-1x\tTotal-PE-aligned-concordantly->1x\tPerc-total-PE-aligned-concordantly->1x\tTotal-PE-aligned-concordantly-or-discordantly-0x\tTotal-PE-aligned-concordantly-or-discordantly-0x-Mates-pair\tTotal-PE-aligned-concordantly-or-discordantly-0x-aligned-0x\tPerc-total-PE-aligned-concordantly-or-discordantly-0x-aligned-0x\tTotal-PE-aligned-concordantly-or-discordantly-0x-aligned-1x\tPerc-total-PE-aligned-concordantly-or-discordantly-0x-aligned-1x\tTotal-PE-aligned-concordantly-or-discordantly-0x-aligned->1x\tPerc-total-PE-aligned-concordantly-or-discordantly-0x-aligned->1x\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

col01="${step}";                                                                             # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f9`;                                         # Sample
col03=`tail -n16 ${pbs_e} | head -n1 | cut -d" " -f1`;                                       # Total-reads
col04=`tail -n2 ${pbs_e} | head -n1 | cut -d" " -f1 | tr -d '%'`;                            # Overal-alignment-rate
col05=`tail -n15 ${pbs_e} | head -n1 | cut -d" " -f3`;                                       # Total-PE
col06=`tail -n15 ${pbs_e} | head -n1 | cut -d" " -f4 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE
col07=`tail -n14 ${pbs_e} | head -n1 | cut -d" " -f5`;                                       # Total-PE-aligned-concordantly-0x
col08=`tail -n14 ${pbs_e} | head -n1 | cut -d" " -f6 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE-aligned-concordantly-0x
col09=`tail -n9 ${pbs_e} | head -n1 | cut -d" " -f7`;                                        # Total-PE-aligned-concordantly-0x-aligned-discordantly-1x
col10=`tail -n9 ${pbs_e} | head -n1 | cut -d" " -f8 | tr -d '(' | tr -d ')' | tr -d '%'`;    # Perc-total-PE-aligned-concordantly-0x-aligned-discordantly-1x
col11=`tail -n13 ${pbs_e} | head -n1 | cut -d" " -f5`;                                       # Total-PE-aligned-concordantly-1x
col12=`tail -n13 ${pbs_e} | head -n1 | cut -d" " -f6 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE-aligned-concordantly-1x
col13=`tail -n12 ${pbs_e} | head -n1 | cut -d" " -f5`;                                       # Total-PE-aligned-concordantly->1x
col14=`tail -n12 ${pbs_e} | head -n1 | cut -d" " -f6 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE-aligned-concordantly->1x
col15=`tail -n7 ${pbs_e} | head -n1 | cut -d" " -f5`;                                        # Total-PE-aligned-concordantly-or-discordantly-0x
col16=`tail -n6 ${pbs_e} | head -n1 | cut -d" " -f7`;                                        # Total-PE-aligned-concordantly-or-discordantly-0x-Mates-pair
col17=`tail -n5 ${pbs_e} | head -n1 | cut -d" " -f9`;                                        # Total-PE-aligned-concordantly-or-discordantly-0x-aligned-0x
col18=`tail -n5 ${pbs_e} | head -n1 | cut -d" " -f10 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE-aligned-concordantly-or-discordantly-0x-aligned-0x
col19=`tail -n4 ${pbs_e} | head -n1 | cut -d" " -f9`;                                        # Total-PE-aligned-concordantly-or-discordantly-0x-aligned-1x
col20=`tail -n4 ${pbs_e} | head -n1 | cut -d" " -f10 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE-aligned-concordantly-or-discordantly-0x-aligned-1x
col21=`tail -n3 ${pbs_e} | head -n1 | cut -d" " -f9`;                                        # Total-PE-aligned-concordantly-or-discordantly-0x-aligned->1x
col22=`tail -n3 ${pbs_e} | head -n1 | cut -d" " -f10 | tr -d '(' | tr -d ')' | tr -d '%'`;   # Perc-total-PE-aligned-concordantly-or-discordantly-0x-aligned->1x
col23=`echo "$(basename ${pbs_e})"`;                                                         # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}\t${col20}\t${col21}\t${col22}\t${col23}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 08S1
#................................................

date ## Collect run summary for step 08S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step08S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tUsable-reads\tData-type\tFraction-reads-failed-to-determine\tFraction-reads-fr-secondstrand-1++,1--,2+-,2-+\tFraction-reads-fr-firststrand-1+-,1-+,2++,2--\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

col01="${step}";                                                                                                                       # Step
col02=`echo "$(basename ${pbs_e})" | cut -d"_" -f10`;                                                                                  # Sample
col03=`tail -n1 ${pbs_e} | cut -d" " -f7`;                                                                                             # Usable-reads
col04=`ls ${folder} | grep "${col02}" | head -n1 | while read file; do tail -n4 ${folder}/${file} | head -n1 | cut -d" " -f3 ; done`;  # Data-type
col05=`ls ${folder} | grep "${col02}" | head -n1 | while read file; do tail -n3 ${folder}/${file} | head -n1 | cut -d" " -f7 ; done`;  # Fraction-reads-failed-to-determine
col06=`ls ${folder} | grep "${col02}" | head -n1 | while read file; do tail -n2 ${folder}/${file} | head -n1 | cut -d" " -f7 ; done`;  # Fraction-reads-fr-secondstrand-1++,1--,2+-,2-+
col07=`ls ${folder} | grep "${col02}" | head -n1 | while read file; do tail -n1 ${folder}/${file} | cut -d" " -f7 ; done`;             # Fraction-reads-fr-firststrand-1+-,1-+,2++,2--
col08=`echo "$(basename ${pbs_e})"`;                                                                                                   # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 08S2
#................................................

date ## Collect run summary for step 08S2 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step08S2" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Count\tSum\tAverage\tMedian\tSmaller\tBigger"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

# Spearman correlation:
max=`grep -v "#" ${folder}/spearman-correlation_*.matrix | head -n1 | tr '\t' '\n' | wc -l`; seq 2 ${max} | while read n; do grep -v "#" ${folder}/spearman-correlation_*.matrix | cut -f${n} | head -n1 | while read sample; do grep "$sample" ${folder}/spearman-correlation_*.matrix | cut -f1,${n} >> ${out_path_step10H1_Bash}/tmp_spearman-samples; done; done

#cut -f2 ${out_path_step10H1_Bash}/tmp_spearman-samples | grep -v "_\|1.0000" | awk 'BEGIN { c = 0; sum = 0; }
#  $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; }
#  END { ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ]; } else { median = ( a[c/2] + a[c/2-1] ) / 2; }
#  OFS="\t"; print c"\t"sum"\t"ave"\t"median"\t"a[0]"\t"a[c-1]; }' >> ${summary_output}

cut -f2- ${folder}/spearman-correlation_*.matrix | grep -v "^#\|_" | tr '\t' '\n' | grep -v "^1.0000" | awk 'BEGIN { c = 0; sum = 0; }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; }
  END { ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ]; } else { median = ( a[c/2] + a[c/2-1] ) / 2; }
  OFS="\t"; print c"\t"sum"\t"ave"\t"median"\t"a[0]"\t"a[c-1]; }' >> ${summary_output}

# Pearson correlation:
max=`grep -v "#" ${folder}/pearson-correlation_*.matrix | head -n1 | tr '\t' '\n' | wc -l`; seq 2 ${max} | while read n; do grep -v "#" ${folder}/pearson-correlation_*.matrix | cut -f${n} | head -n1 | while read sample; do grep "$sample" ${folder}/pearson-correlation_*.matrix | cut -f1,${n} >> ${out_path_step10H1_Bash}/tmp_pearson-samples; done; done

#cut -f2 ${out_path_step10H1_Bash}/tmp_pearson-samples | grep -v "_\|1.0000" | awk 'BEGIN { c = 0; sum = 0; }
#  $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; }
#  END { ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ]; } else { median = ( a[c/2] + a[c/2-1] ) / 2; }
#  OFS="\t"; print c"\t"sum"\t"ave"\t"median"\t"a[0]"\t"a[c-1]; }' >> ${summary_output}

cut -f2- ${folder}/pearson-correlation_*.matrix | grep -v "^#\|_" | tr '\t' '\n' | grep -v "^1.0000" | awk 'BEGIN { c = 0; sum = 0; }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; }
  END { ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ]; } else { median = ( a[c/2] + a[c/2-1] ) / 2; }
  OFS="\t"; print c"\t"sum"\t"ave"\t"median"\t"a[0]"\t"a[c-1]; }' >> ${summary_output}

echo -e "Step\tSample" >> ${out_path_step10H1_Bash}/tmp_info
#echo -e "08S2\tSpearman-replicates" >> ${out_path_step10H1_Bash}/tmp_info
echo -e "${step}\tSpearman" >> ${out_path_step10H1_Bash}/tmp_info
#echo -e "08S2\tPearson-replicates" >> ${out_path_step10H1_Bash}/tmp_info
echo -e "${step}\tPearson" >> ${out_path_step10H1_Bash}/tmp_info

paste ${out_path_step10H1_Bash}/tmp_info ${summary_output} >> ${out_path_step10H1_Bash}/tmp_complete ; mv ${out_path_step10H1_Bash}/tmp_complete ${summary_output}

rm -rf ${out_path_step10H1_Bash}/tmp_info ${out_path_step10H1_Bash}/tmp_pearson-samples ${out_path_step10H1_Bash}/tmp_spearman-samples

fi; done

#................................................
#  Summary 09L1
#................................................

date ## Collect run summary for step 09L1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2`  ; if [[ "${step}" == "step09L1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tn_threshold\tTotal\tLongest\tSmallest\tn_plus\tup_t_n\t0_1k\t1k\t1k_2k\t2k_5k\t5k_10k\t10k_20k\t20k_30k\t30k_40k\t40k_50k\t50kplus\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

col01="${step}";                                        # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f12`;   # Sample
col03=`tail -n16 ${pbs_o} | head -n1 | cut -f2`;        # n_threshold
col04=`tail -n15 ${pbs_o} | head -n1 | cut -f2`;        # Total
col05=`tail -n14 ${pbs_o} | head -n1 | cut -f2`;        # Longest
col06=`tail -n13 ${pbs_o} | head -n1 | cut -f2`;        # Smallest
col07=`tail -n12 ${pbs_o} | head -n1 | cut -f2`;        # n_plus
col08=`tail -n11 ${pbs_o} | head -n1 | cut -f2`;        # up_t_n
col09=`tail -n10 ${pbs_o} | head -n1 | cut -f2`;        # 0_1k
col10=`tail -n09 ${pbs_o} | head -n1 | cut -f2`;        # 1k
col11=`tail -n08 ${pbs_o} | head -n1 | cut -f2`;        # 1k_2k
col12=`tail -n07 ${pbs_o} | head -n1 | cut -f2`;        # 2k_5k
col13=`tail -n06 ${pbs_o} | head -n1 | cut -f2`;        # 5k_10k
col14=`tail -n05 ${pbs_o} | head -n1 | cut -f2`;        # 10k_20k
col15=`tail -n04 ${pbs_o} | head -n1 | cut -f2`;        # 20k_30k
col16=`tail -n03 ${pbs_o} | head -n1 | cut -f2`;        # 30k_40k
col17=`tail -n02 ${pbs_o} | head -n1 | cut -f2`;        # 40k_50k
col18=`tail -n01 ${pbs_o} | cut -f2`;                   # 50kplus
col19=`echo "$(basename ${pbs_o})"`;                    # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}\t${col16}\t${col17}\t${col18}\t${col19}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 09S1
#................................................

date ## Collect run summary for step 09S1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | cut -d"_" -f2` ; n=`echo ${folder} | cut -d"_" -f2 | cut -d"p" -f2` ; if [[ "${step}" == "step09S1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}"` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${n}"` ; out_file=`cut -f4 ${input} | grep "${n}"`
pbs_file=`cut -f5 ${input} | grep "${n}"`
header="Step\tSample\tLane_Read\t#Reads\tAveragePhred\tMedianPhred\tMinPhred\tMaxPhred\tBase-file"
summary_output="${out_path_step10H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

pbs=`echo ${pbs_o} | cut -d"." -f1`
inpfolder=`tail -n1 ${pbs}.pbs | cut -d" " -f5`

sample=`echo "$(basename ${pbs_o})" | cut -d"_" -f9`;

ls ${inpfolder}/*${sample}*zip | cut -d"." -f1-2| while read file; do 
base_file=`echo "$(basename ${file})"`;
unzip -p ${file}.zip ${base_file}/fastqc_data.txt > fastqc_data.txt ;
lines=`grep -n ">>END_MODULE" fastqc_data.txt | cut -d":" -f1 | head -n2 | tail -n1`; y=$((lines-14)) ; z=$((y-1));
grep -A${y} ">>Per base sequence quality" fastqc_data.txt | tail -n${z} > temp.txt

col01="${step}";                                                                                # Step
col02=`echo "$(basename ${pbs_o})" | cut -d"_" -f9`;                                            # Sample
col03=`grep "Filename" fastqc_data.txt | cut -f2 | cut -d"_" -f8 | cut -d"." -f1`;              # Lane_Read
col04=`grep "Total Sequences" fastqc_data.txt | cut -f2`;                                       # #Reads
col05=`cut -f2 temp.txt | awk '{s+=$1} END {print (NR ? s/NR : "NaN")}'`;                       # AveragePhred
col06=`cut -f2 temp.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`;                    # MedianPhred
col07=`cut -f2 temp.txt | awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }'`;  # MinPhred
col08=`cut -f2 temp.txt | awk '{if (NR == 1 || $1 > max) max = $1} END {print max}'`;           # MaxPhred
col09=`echo "$(basename ${pbs_o})"`;                                                            # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}" >> ${summary_output};

done ; rm -rf fastqc_data.txt temp.txt

done; 

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
sed -i 's,${out_path_step10H1_Bash},'"${out_path_step10H1_Bash}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v