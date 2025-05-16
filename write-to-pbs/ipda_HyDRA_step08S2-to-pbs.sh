#!/bin/bash

version="1.0.1"
usage(){
echo "
Written by Isabela Almeida
Created on Jun 08, 2023
Last modified on May 15, 2025
Version: ${version}

Description: Write and submit PBS jobs for Step 08S2 (short-reads) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HyDRA_step08S2-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 1 -c 8 -w "02:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            path/from/working/dir/to/hydra07S1_gen-map_Bowtie2_DATE/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem
                            of coordsorted.bam files and no full stops.

                            Col2:
                            path/from/working/dir/to/hydra07S1_gen-map_Bowtie2_DATE/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_coordsorted.bam

                            Col3:
                            Semicolon separated hex colors or name colors - one of each sample in your dataset.
                            Note: Same string should be provided for all lines in the input for a single run.
                            e.g. for n samples 30, Col3 containg either
                            #61bcc6;#ade66c;#f11e71;#9cd5d8;#6213fe;#e17ebb;#f7dc02;#a14b74;#a66fc7;#de08e3;#be1f5f;#1838f9;#d6fd2e;#16a184;#5b1119;#b9f9f7;#57e148;#0c4b0e;#a88f99;#4b8623;#0c170d;#cfc154;#424cab;#669dd7;#3a06a8;#ec5a5e;#56dac1;#e4d2d6;#262f4f;#b3d17e
                            OR
                            black;grey;lightgrey;palevioletred;firebrick;crimson;orangered;saddlebrown;bisque;navajowhite;orchid;magenta;darkkhaki;olive;darkolivegreen;darkseagreen;darkgreen;mediumseagreen;aquamarine;lightcyan;darkcyan;cadetblue;lightskyblue;lightslategrey;royalblue;darkblue;darkslateblue;indigo;plum;fuchsia;hotpink;lightpink

                            Col4:
                            Space separated labels for each sample - same order as files provided in Col2 of entire file
                            Note: Same string should be provided for all lines in the input for a single run.

                            It does not matter if same stem 
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
#-->08 [S] Assess strandness (1RSeQC) & Assess replicates correlation (2DeepTools)
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
        i) input="${OPTARG}";;       # Input files including path
        p) pbs_stem="${OPTARG}";;    # Stem for PBS file names
        e) email="${OPTARG}";;       # Email for PBS job
        m) mem="${OPTARG}";;         # Memory required for PBS job
        c) ncpus="${OPTARG}";;       # Number of CPUS required for PBS job
        w) walltime="${OPTARG}";;    # Clock walltime required for PBS job
        h) Help ; exit;;             # Print Help and exit
        v) echo "${version}"; exit;; # Print version and exit
        ?) echo script usage: bash ipda_HyDRA_step08S2-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_hydra08S2-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# Python3
# <https://deeptools.readthedocs.io/en/develop/content/installation.html>
module_deeptools=conda-envs/deeptools-3.5.4

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step08S2_DeepTools="hydra08S2_rep-cor_DeepTools_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step08S2_DeepTools}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HyDRA_step08S2-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step08S2_DeepTools}"
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
echo "## Executing bash ipda_HyDRA_step08S2-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step08S2_DeepTools}"
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
echo "module load ${module_deeptools}" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write PBS command lines
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Run step" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "## Create color palette at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
hex_name_colors=`cut -f3 ${input} | sort | uniq`; echo "echo '${hex_name_colors}' | tr ';' '\n' >> ${out_path_step08S2_DeepTools}/pca-color-pallete_${thislogdate}.txt" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Create PCA labels at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
labels=`cut -f4 ${input} | sort | uniq`; echo "echo '${labels}' | tr ' ' '\n' | uniq >> ${out_path_step08S2_DeepTools}/pca-labels_${thislogdate}.txt" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Create color scheme for PCA plot at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
labels=`cut -f4 ${input} | sort | uniq`; sample_num=`echo ${labels} | tr ' ' '\n' | uniq | awk 'END{print NR}'` ; echo "head -n${sample_num} ${out_path_step08S2_DeepTools}/pca-color-pallete_${thislogdate}.txt >> ${out_path_step08S2_DeepTools}/pca-color-scheme_${thislogdate}.txt" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Create label-color scheme for PCA plot at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "paste ${out_path_step08S2_DeepTools}/pca-labels_${thislogdate}.txt ${out_path_step08S2_DeepTools}/pca-color-scheme_${thislogdate}.txt >> ${out_path_step08S2_DeepTools}/pca-labels-color-scheme_${thislogdate}.txt" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Separate colors from label-color scheme at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "cut -f2 ${out_path_step08S2_DeepTools}/pca-labels-color-scheme_${thislogdate}.txt >> ${out_path_step08S2_DeepTools}/pca-sample-colors_${thislogdate}.txt" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Generate BAM comparison summary | Run multiBamSummary at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
bams=`cut -f2 ${input} | sort | uniq | tr '\n' ' '`; labels=`cut -f4 ${input} | sort | uniq`; echo "multiBamSummary bins --bamfiles ${bams} --labels ${labels} --outFileName ${out_path_step08S2_DeepTools}/bam-comparison-summary_${pbs_stem}.npz --maxFragmentLength 280 --numberOfProcessors ${ncpus}" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Generate Spearman correlation | Run plotCorrelation at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
labels=`cut -f4 ${input} | sort | uniq`; echo "plotCorrelation --corData ${out_path_step08S2_DeepTools}/bam-comparison-summary_${pbs_stem}.npz --corMethod spearman --whatToPlot heatmap --labels ${labels} --plotFile ${out_path_step08S2_DeepTools}/spearman-correlation_${pbs_stem}.pdf --plotTitle Spearman-Correlation_${pbs_stem}_${thislogdate} --outFileCorMatrix ${out_path_step08S2_DeepTools}/spearman-correlation_${pbs_stem}.matrix --plotNumbers" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Generate Pearson correlation | Run plotCorrelation at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
labels=`cut -f4 ${input} | sort | uniq`; echo "plotCorrelation --corData ${out_path_step08S2_DeepTools}/bam-comparison-summary_${pbs_stem}.npz --corMethod pearson --whatToPlot heatmap --labels ${labels} --plotFile ${out_path_step08S2_DeepTools}/pearson-correlation_${pbs_stem}.pdf --plotTitle Pearson-Correlation_${pbs_stem}_${thislogdate} --outFileCorMatrix ${out_path_step08S2_DeepTools}/pearson-correlation_${pbs_stem}.matrix --plotNumbers" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

echo 'echo "## Generate PCA plot | Run plotPCA at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
labels=`cut -f4 ${input} | sort | uniq`; sample_num=`echo ${labels} | tr ' ' '\n' | uniq | awk 'END{print NR}'` ; hex_name_colors=`cut -f3 ${input} | sort | uniq`; colors=`echo ${hex_name_colors} | tr ';' '\n' | tr '\n' ' ' | head -n${sample_num}` ; echo "plotPCA --corData ${out_path_step08S2_DeepTools}/bam-comparison-summary_${pbs_stem}.npz --plotHeight 20 --labels ${labels} --colors ${colors} --plotFile ${out_path_step08S2_DeepTools}/pca-plot_${pbs_stem}.pdf --plotTitle PCA-plot_${pbs_stem}_${thislogdate} --log2" >> ${pbs_stem}_${thislogdate}.pbs

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n25 ${pbs} | head -n1)" ; echo "#                       $(tail -n22 ${pbs} | head -n1)" ; echo "#                       $(tail -n19 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n13 ${pbs} | head -n1)" ; echo "#                       $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 08S2 jobs) at
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
sed -i 's,${module_deeptools},'"${module_deeptools}"',g' "$logfile"
sed -i 's,${out_path_step08S2_DeepTools},'"${out_path_step08S2_DeepTools}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v