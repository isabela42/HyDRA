#!/bin/bash

version="1.0.1"
usage(){
echo "
Written by Isabela Almeida
Created on Feb 11, 2024
Last modified on May 15, 2025
Version: ${version}

Description: Write and submit PBS jobs for Step 17H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HyDRA_step17H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "-9.email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 1 -c 1 -w "24:00:00"

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
#-->17 [H] Summary of hybrid de novo transcriptome assembly steps

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
        ?) echo script usage: bash ipda_HyDRA_step17H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_hydra17H1-to-pbs_${thislogdate}.txt

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step17H1_Bash="hydra17H1_summary_Bash_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step17H1_Bash}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HyDRA_step17H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step17H1_Bash}"
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
echo "## Executing bash ipda_HyDRA_step17H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step17H1_Bash}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Summary 12H1
#................................................

date ## Collect run summary for step 12H1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | rev | cut -d'/' -f1 | rev | cut -d'_' -f1` ; if [[ "${step}" == "hydra12H1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}" | rev | cut -d'/' -f1 | rev` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${step}"` ; out_file=`cut -f4 ${input} | grep "${step}"`
pbs_file=`cut -f5 ${input} | grep "${step}"`
header="Step\tSample\tComplete\tComplete-perc\tComplete-sc\tComplete-sc-perc\tComplete-duplicated\tComplete-duplicated-perc\tFragmented\tFragmented-perc\tMissing\tMissing-perc\tTotal\tBase-file"
summary_output="${out_path_step17H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

col01="${step}";                                                                           # Step
col02=`echo "$(basename ${pbs_o})" | rev | cut -d"_" -f2- | rev | cut -d'_' -f4-`;         # Sample
col03=`grep -A2 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d" " -f1`;                 # Complete
col04=`grep -A1 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d"%" -f1 | cut -d":" -f2`; # Complete-perc
col05=`grep -A3 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d" " -f1`;                 # Complete-sc
col06=`grep -A1 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d"%" -f2 | cut -d":" -f2`; # Complete-sc-perc
col07=`grep -A4 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d" " -f1`;                 # Complete-duplicated
col08=`grep -A1 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d"%" -f3 | cut -d":" -f2`; # Complete-duplicated-perc
col09=`grep -A5 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d" " -f1`;                 # Fragmented
col10=`grep -A1 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d"%" -f4 | cut -d":" -f2`; # Fragmented-perc
col11=`grep -A6 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d" " -f1`;                 # Missing
col12=`grep -A1 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d"%" -f5 | cut -d":" -f2`; # Missing-perc
col13=`grep -A7 "Results:" ${pbs_o} | tail -n1 | cut -f2 | cut -d" " -f1`;                 # Total
col14=`echo "$(basename ${pbs_o})"`;                                                       # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 12H2
#................................................

date ## Collect run summary for step 12H2 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | rev | cut -d'/' -f1 | rev | cut -d'_' -f1` ; if [[ "${step}" == "hydra12H2" ]]; then

run_name=`cut -f1 ${input} | grep "${step}" | rev | cut -d'/' -f1 | rev` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${step}"` ; out_file=`cut -f4 ${input} | grep "${step}"`
pbs_file=`cut -f5 ${input} | grep "${step}"`
header="Step\tSample\t#Genes\t#Transcripts\tGC%\tN10\tN20\tN30\tN40\tN50\tMedian-len\tAverage-len\tAssembled-bases\tBase-file"
summary_output="${out_path_step17H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${out_file} | while read pbs_o ; do

col01="${step}";                                                                                                                         # Step
col02=`echo "$(basename ${pbs_o})" | rev | cut -d"_" -f2- | rev | cut -d'_' -f4-`;                                                       # Sample
col03=`ls ${folder} | grep "${col02}" | while read file ; do grep "Total trinity 'genes':" ${folder}/${file} | cut -f2 ; done`;          # #Genes
col04=`ls ${folder} | grep "${col02}" | while read file ; do grep "Total trinity transcripts:" ${folder}/${file} | cut -f2 ; done`;      # #Transcripts
col05=`ls ${folder} | grep "${col02}" | while read file ; do grep "Percent GC:" ${folder}/${file} | cut -d' ' -f3 ; done`;                     # GC%
col06=`ls ${folder} | grep "${col02}" | while read file ; do grep -A3 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f3 ; done`;        # N10
col07=`ls ${folder} | grep "${col02}" | while read file ; do grep -A4 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f3 ; done`;        # N20
col08=`ls ${folder} | grep "${col02}" | while read file ; do grep -A5 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f3 ; done`;        # N30
col09=`ls ${folder} | grep "${col02}" | while read file ; do grep -A6 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f3 ; done`;        # N40
col10=`ls ${folder} | grep "${col02}" | while read file ; do grep -A7 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f3 ; done`;        # N50
col11=`ls ${folder} | grep "${col02}" | while read file ; do grep -A9 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f4 ; done`;        # Median-len
col12=`ls ${folder} | grep "${col02}" | while read file ; do grep -A10 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f3 ; done`;       # Average-len
col13=`ls ${folder} | grep "${col02}" | while read file ; do grep -A11 "ALL" ${folder}/${file} | tail -n1 | cut -d" " -f4 ; done`;       # Assembled-bases
col14=`echo "$(basename ${pbs_o})"`;                    # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}" >> ${summary_output}; done

fi; done

#................................................
#  Summary 14H1
#................................................

date ## Collect run summary for step 14H1 at
cut -f1 ${input} | while read folder; do
step=`echo ${folder} | rev | cut -d'/' -f1 | rev | cut -d'_' -f1` ; if [[ "${step}" == "hydra14H1" ]]; then

run_name=`cut -f1 ${input} | grep "${step}" | rev | cut -d'/' -f1 | rev` ; log_file=`cut -f2 ${input} | grep "${step}"`
err_file=`cut -f3 ${input} | grep "${step}"` ; out_file=`cut -f4 ${input} | grep "${step}"`
pbs_file=`cut -f5 ${input} | grep "${step}"`
header="Step\tSample\tduplicates-unspliced\tduplicates-spliced\tduplicates-unaligned\tsingle-unspliced\tsingle-spliced\tsingle-unaligned\tNo-path-found\tMap-to-1-loc\tMap-to-2-loc\tMap-to-3-loc\tMap-to-5-10-loc\tMap-to-10+-loc\tBase-file"
summary_output="${out_path_step17H1_Bash}/summary_${run_name}.tsv"
echo -e "${header}" > ${summary_output}

ls ${err_file} | while read pbs_e ; do

col01="${step}";                                                                                                                                                                                                                                         # Step
col02=`echo "$(basename ${pbs_e})" | rev | cut -d"_" -f2- | rev | cut -d'_' -f4-`;                                                                                                                                                                       # Sample #-f5-6 for gems, -f5 for hydra development
col03=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | awk '$1 ~ /@/ || ($6 !~ /N/ && $6 !~ /\*/)' | cut -f1 | sort | uniq -c | cut -d" " -f7- | grep -v "^1" | wc -l ; done`;  # duplicates-unspliced
col04=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | awk '$1 ~ /@/ || ($6 ~ /N/)' | cut -f1 | sort | uniq -c | cut -d" " -f7- | grep -v "^1" | wc -l ; done `;                # duplicates-spliced
col05=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | awk '$1 ~ /@/ || ($6 ~ /\*/)' | cut -f1 | sort | uniq -c | cut -d" " -f7- | grep -v "^1" | wc -l ; done`;                # duplicates-unaligned
col06=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | awk '$1 ~ /@/ || ($6 !~ /N/ && $6 !~ /\*/)' | cut -f1 | sort | uniq -c | cut -d" " -f7- | grep "^1" | wc -l ; done `;    # single-unspliced
col07=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | awk '$1 ~ /@/ || ($6 ~ /N/)' | cut -f1 | sort | uniq -c | cut -d" " -f7- | grep "^1" | wc -l ; done`;                    # single-spliced
col08=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | awk '$1 ~ /@/ || ($6 ~ /\*/)' | cut -f1 | sort | uniq -c | cut -d" " -f7- | grep "^1" | wc -l ; done`;                   # single-unaligned
col09=`grep -c "No paths found for" ${pbs_e}`;                                                                                                                                                                                                           # No-path-found
col10=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 | tr -s " " | tr ' ' '\t' | grep -P "[0-9]\t40" | cut -f2 ; done`; # Map-to-1-loc
col11=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 | tr -s " " | tr ' ' '\t' | grep -P "[0-9]\t3" | cut -f2 ; done`;  # Map-to-2-loc
col12=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 | tr -s " " | tr ' ' '\t' | grep -P "[0-9]\t2" | cut -f2 ; done`;  # Map-to-3-loc
col13=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 | tr -s " " | tr ' ' '\t' | grep -P "[0-9]\t1" | cut -f2 ; done`;  # Map-to-5-10-loc
col14=`ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | cut -f5 | sort | uniq -c | sort -k1,1 -nbr | head -n5 | tr -s " " | tr ' ' '\t' | grep -P "[0-9]\t0" | cut -f2 ; done`;  # Map-to-10+-loc
col15=`echo "$(basename ${pbs_e})"`;                                                                                                                                                                                                                     # Base-file

echo -e "${col01}\t${col02}\t${col03}\t${col04}\t${col05}\t${col06}\t${col07}\t${col08}\t${col09}\t${col10}\t${col11}\t${col12}\t${col13}\t${col14}\t${col15}" >> ${summary_output}

echo -e "Count\tSum\tAverage\tMedian\tSmaller\tBigger" >> ${out_path_step17H1_Bash}/summary_${run_name}_statistics.tsv
ls ${folder} | grep "${col02}" | while read file ; do grep -v "^@" ${folder}/${col02}/*_splicing-assessment.genome.sam | cut -f5 | awk 'BEGIN { c = 0; sum = 0; }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ { a[c++] = $1; sum += $1; }
  END { ave = sum / c; if( (c % 2) == 1 ) { median = a[ int(c/2) ]; } else { median = ( a[c/2] + a[c/2-1] ) / 2; }
  OFS="\t"; print c"\t"sum"\t"ave"\t"median"\t"a[0]"\t"a[c-1]; }' >> ${out_path_step17H1_Bash}/summary_${run_name}_statistics.tsv ; done ; done

echo -e "Step\tSample" >> ${out_path_step17H1_Bash}/tmp_info
echo -e "${col01}\t${col02}" >> ${out_path_step17H1_Bash}/tmp_info
paste ${out_path_step17H1_Bash}/summary_${run_name}_statistics.tsv >> ${out_path_step17H1_Bash}/tmp_complete ; mv ${out_path_step17H1_Bash}/tmp_complete ${out_path_step17H1_Bash}/summary_${run_name}_statistics.tsv ; rm -rf ${out_path_step17H1_Bash}/tmp_info

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
sed -i 's,${out_path_step17H1_Bash},'"${out_path_step17H1_Bash}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v