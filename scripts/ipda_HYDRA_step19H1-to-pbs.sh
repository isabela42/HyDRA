#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Jan 08, 2023
Last modified on Jun 23, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 19H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step19H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 120 -c 16 -w "60:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/HYDRA_step15H*_read-support_Samtools_DATE/folder/stem.fasta
                            of .spliced.fasta and .unspliced.fasta files (with full file name)                            

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
#-->19 [H] Identify long noncoding RNAs (1FEELnc)
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
        ?) echo script usage: bash ipda_HYDRA_step19H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step19H1-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

# Genome file from Gencode.
# Donwloaded from:
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz
gencodegenome=/working/lab_julietF/isabelaA/data/references/GRCh38.p13.genome.fa

# Annotation file for all human transcripts from Gencode.
# BED format, generated above.
gencodebed=/working/lab_julietF/isabelaA/data/references/gencode.v36.annotation.bed

# Annotation files for human transcriptome.
# Genome build GRCh38.p2 updated in 2015.
# Contains 65309 annotated genes (being 92 ERCC entries)
transgtf=/working/lab_julietF/isabelaA/data/references/TranscriptomeGRCh38rel79_ERCC.gtf
#transplusgtf=/working/lab_julietF/isabelaA/data/references/TranscriptomeGRCh38rel79_plus_mencRNAs.gtf
transcriptome=/working/lab_julietF/isabelaA/data/references/TranscriptomeGRCh38p7_ERCC.fa

# Annotation file for human protein-coding transcripts from Gencode.
# Obtained from above, retrieving only "protein_coding" type entries
# Protein-coding transcripts were retrieved using grep "^#\|protein_coding"
ptngtf=/working/lab_julietF/isabelaA/data/references/gencode.v36.proteincoding.gtf

# Annotation file for human lncRNAs from Gencode.
# Downloaded from:
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.long_noncoding_RNAs.gtf.gz
lncrnagtf=/working/lab_julietF/isabelaA/data/references/gencode.v36.long_noncoding_RNAs.gtf
# grep -v "TEC" out of above (to be experimentally confirmed transcripts)
goodlncrnagtf=/working/lab_julietF/isabelaA/data/references/gencode.v36.confirmed_long_noncoding_RNAs.gtf

#................................................
#  Required modules, softwares and libraries
#................................................

## Load tools from HPC
# For more info, see
# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>

# UCSC Tools 20160223 (genePredToGtf and bedToGenePred)
module_ucsctools=ucsctools/20160223

# BedTools 2.29.0
module_bedtools=bedtools/2.29.0

# Samtools 1.9
module_samtools=samtools/1.9

# Minimap2 2.16
module_minimap2=minimap2/2.16

# R 3.5.1
module_R=R/3.5.1

# FEELnc
module_feelnc=FEELnc/0.2
feelnc_path=/software/FEELnc/FEELnc-v.0.2 # from which feelnc

## Path to user-installed tools

# NA

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step19H1_FEELnc="HYDRA_step19H1_lncRNA-potential_FEELnc_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step19H1_FEELnc}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step19H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step19H1_FEELnc}"
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
echo "## Executing bash ipda_HYDRA_step19H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step19H1_FEELnc}"
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
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Load tools from HPC"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_ucsctools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_bedtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_samtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_minimap2}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_R}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_feelnc}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "export FEELNCPATH=${feelnc_path}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'export PATH="'"${feelnc_path}"'/FEELnc/bin/LINUX:$PATH"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Path to user-installed tools"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "# None required"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Align transcripts to genome with Minimap2 at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; mkdir -p ${out_path_step19H1_FEELnc}/${file} ; echo "minimap2 -x splice --secondary=no -a -t ${ncpus} -o ${out_path_step19H1_FEELnc}/${file}/${file}.sam ${gencodegenome} ${path_file}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Sort alignment files and save as bam at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "samtools view -u ${out_path_step19H1_FEELnc}/${file}/${file}.sam | samtools sort -@ ${ncpus} -o ${out_path_step19H1_FEELnc}/${file}/${file}.sorted.bam" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Transform file to bed format at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bamToBed -splitD -bed12 -i ${out_path_step19H1_FEELnc}/${file}/${file}.sorted.bam > ${out_path_step19H1_FEELnc}/${file}/${file}.bed12" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Intersect alignments with Gencode annotation for future use at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -wo -s -f 0.75 -r -a ${out_path_step19H1_FEELnc}/${file}/${file}.bed12 -b ${gencodebed} >> ${out_path_step19H1_FEELnc}/${file}/${file}.gencode.genomefeatures.bed " >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Create gtf file for FEELnc at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedToGenePred ${out_path_step19H1_FEELnc}/${file}/${file}.bed12 stdout | genePredToGtf file stdin ${out_path_step19H1_FEELnc}/${file}/${file}.gtf " >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run FEELnc Step 1 - FEELnc_filter.pl at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "FEELnc_filter.pl -i ${out_path_step19H1_FEELnc}/${file}/${file}.gtf -a $transgtf --biotype transcript_biotype=protein_coding --size=200 --monoex=0 -p ${ncpus} -o ${out_path_step19H1_FEELnc}/${file}/${file}.feelncfilter.log > ${out_path_step19H1_FEELnc}/${file}/${file}.lncRNAcandidates.gtf" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
# Additional parameters to consider: --minfrac_over=X | minimal fraction out of the candidate lncRNA size to be considered for overlap (default = 0)

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run FEELnc Step 2 - FEELnc_codpot.pl at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "FEELnc_codpot.pl -i ${out_path_step19H1_FEELnc}/${file}/${file}.lncRNAcandidates.gtf -a $ptngtf -l $goodlncrnagtf -g $gencodegenome -o ${file} --outdir=${out_path_step19H1_FEELnc}/${file}/codpot_${file}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run FEELnc Step 3 - FEELnc_classifier.pl at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "FEELnc_classifier.pl -i ${out_path_step19H1_FEELnc}/${file}/codpot_${file}/${file}.lncRNA.gtf -a $ptngtf -l ${out_path_step19H1_FEELnc}/${file}/codpot_${file}/${file}.classifier.log >> ${out_path_step19H1_FEELnc}/${file}/codpot_${file}/${file}.lncRNAclasses.txt" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Move figures at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; mkdir -p ${out_path_step19H1_FEELnc}/VarImpPlots ${out_path_step19H1_FEELnc}/ROCcurves ; echo "mv ${out_path_step19H1_FEELnc}/${file}/codpot_${file}/*varImpPlot.png ${out_path_step19H1_FEELnc}/VarImpPlots ; mv ${out_path_step19H1_FEELnc}/${file}/codpot_${file}/*TGROC.png ${out_path_step19H1_FEELnc}/ROCcurves" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n28 ${pbs} | head -n1)" ; echo "#                       $(tail -n25 ${pbs} | head -n1)" ; echo "#                       $(tail -n22 ${pbs} | head -n1)" ; echo "#                       $(tail -n19 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n13 ${pbs} | head -n1)" ; echo "#                       $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done


date ## Status of all user jobs (including HYDRA step 19H1 jobs) at
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
sed -i 's,${module_ucsctools},'"${module_ucsctools}"',g' "$logfile"
sed -i 's,${module_bedtools},'"${module_bedtools}"',g' "$logfile"
sed -i 's,${module_samtools},'"${module_samtools}"',g' "$logfile"
sed -i 's,${module_minimap2},'"${module_minimap2}"',g' "$logfile"
sed -i 's,${module_R},'"${module_R}"',g' "$logfile"
sed -i 's,${module_feelnc},'"${module_feelnc}"',g' "$logfile"
sed -i 's,${feelnc_path},'"${feelnc_path}"',g' "$logfile"
sed -i 's,${gencodegenome},'"${gencodegenome}"',g' "$logfile"
sed -i 's,${gencodebed},'"${gencodebed}"',g' "$logfile"
sed -i 's,${transgtf},'"${transgtf}"',g' "$logfile"
sed -i 's,${ptngtf},'"${ptngtf}"',g' "$logfile"
sed -i 's,${lncrnagtf},'"${lncrnagtf}"',g' "$logfile"
sed -i 's,${goodlncrnagtf},'"${goodlncrnagtf}"',g' "$logfile"
sed -i 's,${out_path_step19H1_FEELnc},'"${out_path_step19H1_FEELnc}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v