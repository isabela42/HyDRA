#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Oct 22, 2023
Last modified on Jun 23, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 16H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step16H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 20 -c 16 -w "10:00:00"

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
#-->16 [H] Annotation (1GMAP, 1sam2bed, 1BedTools)
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
        ?) echo script usage: bash ipda_HYDRA_step16H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step16H1-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

# GEMS
reference_gmap_db=GenomeGRCh38p7_gmap
referencedir=/working/lab_julietF/isabelaA/data/references/
stem_genome="GenomeGRCh38p7"

# Transcriptome BED file:
transcriptome=/working/lab_julietF/mainaB/ReferenceGenomes/TranscriptID_TranscriptomeGRCh38rel79.bed

# In-House database of lncRNA genes (previously generated):
lncrnadb=/working/lab_julietF/isabelaA/data/references/InHouse_lncRNAdb.bed


# # HYDRA development
# # Reference genome/transcriptome:
# reference_gmap_db=GenomeT2TCHM13v2_gmap_db # available at /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/GenomeT2TCHM13v2_gmap_db
# referencedir=/working/lab_julietF/isabelaA/data/references/T2TCHM13v2/
# stem_genome="T2TCHM13v2"

# # Transcriptome BED file:
# # This file was created from
# # /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes.bed
# # pg parent gene according to t2t ; gid gene id from t2t ; gbt gene biotype from tt; gg gene id of pg in hg38 ; g  gene biotype of pg in hg38; if hg38 not found info becomes the one from T2T
# # with the following command:
# # grep -P "ensembl\ttranscript" /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes.bed | while read line; do pt=`echo ${line} | grep -Po 'projection_parent_transcript "\K.*?(?=")' | cut -d"." -f1` ; pg=`echo ${line} | grep -Po 'projection_parent_gene "\K.*?(?=")' | cut -d"." -f1` ; tid=`echo ${line} | grep -Po 'transcript_id "\K.*?(?=")'` ; gid=`echo ${line} | grep -Po 'gene_id "\K.*?(?=")'` ; gbt=`echo ${line} | grep -Po 'gene_biotype "\K.*?(?=")'` ; tbt=`echo ${line} | grep -Po 'transcript_biotype "\K.*?(?=")'` ; if grep -q "gene_id \"${pg}\".*transcript_id \"${pt}\"" /working/lab_julietF/mainaB/ReferenceGenomes/TranscriptID_TranscriptomeGRCh38rel79.bed ; then grep "gene_id \"${pg}\".*transcript_id \"${pt}\"" /working/lab_julietF/mainaB/ReferenceGenomes/TranscriptID_TranscriptomeGRCh38rel79.bed | while read found; do g=`echo ${found} | grep -Po 'gene_biotype "\K.*?(?=")'` ; writeg=$g ; t=`echo ${found} | grep -Po 'transcript_biotype "\K.*?(?=")'` ; writet=$t ; gg=`echo ${found} | grep -Po 'gene_id "\K.*?(?=")'` ; infogg="hg38gene" ; writegg=$gg ; tt=`echo ${found} | grep -Po 'transcript_id "\K.*?(?=")'` ; infott="hg38transcript" ; writett=$tt ; echo -e "${line}\tT2Tgene\t${gid}\t${gbt}\t${infogg}\t${writegg}\t${writeg}\tT2Ttranscript\t${tid}\t${tbt}\t${infott}\t${writett}\t${writet}" >> /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_hg38annotation_gt.bed ; done ; else infogg="hg38geneidnotfound" ; infott="hg38transcriptidnotfound" ; writegg=$gid ; writett=$tid ; writeg=$gbt ; writet=$tbt ; echo -e "${line}\tT2Tgene\t${gid}\t${gbt}\t${infogg}\t${writegg}\t${writeg}\tT2Ttranscript\t${tid}\t${tbt}\t${infott}\t${writett}\t${writet}" >> /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_hg38annotation_gt.bed ; fi ; done
# # grep -P "ensembl\ttranscript" /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes.bed | while read line; do pt=`echo ${line} | grep -Po 'projection_parent_transcript "\K.*?(?=")' | cut -d"." -f1` ; pg=`echo ${line} | grep -Po 'projection_parent_gene "\K.*?(?=")' | cut -d"." -f1` ; tid=`echo ${line} | grep -Po 'transcript_id "\K.*?(?=")'` ; gid=`echo ${line} | grep -Po 'gene_id "\K.*?(?=")'` ; gbt=`echo ${line} | grep -Po 'gene_biotype "\K.*?(?=")'` ; tbt=`echo ${line} | grep -Po 'transcript_biotype "\K.*?(?=")'` ; if grep -q "transcript_id \"${pt}\"" /working/lab_julietF/mainaB/ReferenceGenomes/TranscriptID_TranscriptomeGRCh38rel79.bed ; then grep "transcript_id \"${pt}\"" /working/lab_julietF/mainaB/ReferenceGenomes/TranscriptID_TranscriptomeGRCh38rel79.bed | while read found; do g=`echo ${found} | grep -Po 'gene_biotype "\K.*?(?=")'` ; writeg=$g ; t=`echo ${found} | grep -Po 'transcript_biotype "\K.*?(?=")'` ; writet=$t ; gg=`echo ${found} | grep -Po 'gene_id "\K.*?(?=")'` ; infogg="hg38gene" ; writegg=$gg ; tt=`echo ${found} | grep -Po 'transcript_id "\K.*?(?=")'` ; infott="hg38transcript" ; writett=$tt ; echo -e "${line}\tT2Tgene\t${gid}\t${gbt}\t${infogg}\t${writegg}\t${writeg}\tT2Ttranscript\t${tid}\t${tbt}\t${infott}\t${writett}\t${writet}" >> /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_hg38annotation_t.bed ; done ; else infogg="hg38geneidnotfound" ; infott="hg38transcriptidnotfound" ; writegg=$gid ; writett=$tid ; writeg=$gbt ; writet=$tbt ; echo -e "${line}\tT2Tgene\t${gid}\t${gbt}\t${infogg}\t${writegg}\t${writeg}\tT2Ttranscript\t${tid}\t${tbt}\t${infott}\t${writett}\t${writet}" >> /working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_hg38annotation_t.bed ; fi ; done
# transcriptome=/working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-genes_hg38annotation_t.bed

# # In-House database of lncRNA genes (previously generated):
# lncrnadb=/working/lab_julietF/isabelaA/data/references/InHouse_lncRNAdb.bed

#................................................
#  Required modules, softwares and libraries
#................................................

## Load tools from HPC
# For more info, see
# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>

# Bedops 2.4.41
# Larger inputs only run with megarow: <https://www.biostars.org/p/283321/>
# Larger inputs also require not sorting: <https://github.com/bedops/bedops/issues/208>
module_bedops=bedops/2.4.41

# SamTools 1.9
module_samtools=samtools/1.9

# BedTools 2.29.0
module_bedtools=bedtools/2.29.0

# GMAP 2023-07-20
module_gmap=gmap/20230720

## Path to user-installed tools

# None required

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step16H1_BedTools="HYDRA_step16H1_transcriptome-annotation_BedTools_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step16H1_BedTools}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step16H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step16H1_BedTools}"
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
echo "## Executing bash ipda_HYDRA_step16H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step16H1_BedTools}"
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
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_bedops}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_samtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_bedtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "module load ${module_gmap}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Path to user-installed tools"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "# None required"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Run GMAP alignment transcripts to genome at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; mkdir -p ${out_path_step16H1_BedTools}/${file} ; echo "gmap -n 1 -t ${ncpus} -B 5 --dir=${referencedir} --db=${reference_gmap_db} ${path_file}  --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse > ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.sam" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Convert SAM to BED at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "convert2bed-megarow --input=sam --do-not-sort --max-mem=${mem}G < ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.sam > ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.bed12" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Get main BED columns at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f1-6 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.bed12 >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Get genome features with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -wo -s -f 0.75 -r -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.bed -b ${transcriptome} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Get genome counts with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -s -f 0.75 -r -C -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.bed -b ${transcriptome} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomecounts.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Identify transcripts that are not in the genome with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -v -wa -s -f 0.75 -r -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.bed -b ${transcriptome} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.notingenome.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Get lncRNA DB features with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -wo -s -f 0.75 -r -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.notingenome.bed -b ${lncrnadb} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.lncRNAfeatures.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Get lncRNA DB counts with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -s -f 0.75 -r -C -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.notingenome.bed -b ${lncrnadb} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.lncRNAcounts.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Identify novel transcripts with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect -v -wa -s -f 0.75 -r -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.notingenome.bed -b ${lncrnadb} ${transcriptome} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.novel.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Identify potential fragment on novel transcripts with BedTools at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "bedtools intersect  -wo -s -f 0.5 -e -a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.novel.bed -b ${lncrnadb} ${transcriptome} >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.potentialfragment.bed" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## Create annotation summary at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
# Next line is for a normal ensembl annotation file
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f4 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_a ; cut -d'\"' -f2,6 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed | tr '\"' '\t' >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_b ; grep -Po 'gene_biotype \"\K.*?(?=\")' ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_c ; grep -Po 'transcript_biotype \"\K.*?(?=\")' ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_d ; paste ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_b ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_c ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_d >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation; rm -rf ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_b ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_c ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_d ; cut -f4,10,11 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.lncRNAfeatures.bed | sed \"s/\t./\.\tncRNAdb\t\./g\" >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
# Next line is for the T2T and hg38 combine annotation file
# cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f4,18,24 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_a ; cut -f22,28 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_b ; paste ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_b >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation; rm -rf ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_a ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.tmp_b ; cut -f4,10,11 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.lncRNAfeatures.bed | sed \"s/\t./\.\tncRNAdb\t\./g\" >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f1,4 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation | sort | uniq >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "grep \"protein_coding\" ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary.protein.transcriptID" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "grep \"lincRNA\|ncRNAdb\|lncRNA\" ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary.ncRNA.transcriptID" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "grep \"antisense\" ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary.antisense.transcriptID" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "grep \"pseudogene\" ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.annotation.summary.pseudogene.transcriptID" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo 'echo "## List ID and type of retrieved annotated genes at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -f12,16 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed | tr '\"' '\t' | sort | uniq >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}_geneID-type_retrieved-annotated-genes" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
# For common ensembl annotation
#cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file}" | sed 's/\(.*\)\..*/\1/')"` ; echo "cut -d'\"' -f2,14 ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}.genomefeatures.bed | tr '\"' '\t' | sort | uniq >> ${out_path_step16H1_BedTools}/${file}/${file}.${stem_genome}_geneID-type_retrieved-annotated-genes" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n39 ${pbs} | head -n1)" ; echo "#                       $(tail -n36 ${pbs} | head -n1)" ;  echo "#                       $(tail -n33 ${pbs} | head -n1)" ; echo "#                       $(tail -n30 ${pbs} | head -n1)" ; echo "#                       $(tail -n27 ${pbs} | head -n1)" ; echo "#                       $(tail -n24 ${pbs} | head -n1)" ; echo "#                       $(tail -n21 ${pbs} | head -n1)" ; echo "#                       $(tail -n18 ${pbs} | head -n1)" ; echo "#                       $(tail -n15 ${pbs} | head -n1)" ;  echo "#                       $(tail -n12 ${pbs} | head -n1)" ;  echo "#                       $(tail -n9 ${pbs} | head -n1)" ;  echo "#                       $(tail -n8 ${pbs} | head -n1)" ;  echo "#                       $(tail -n7 ${pbs} | head -n1)" ;  echo "#                       $(tail -n6 ${pbs} | head -n1)" ; echo "#                       $(tail -n5 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 16H1 jobs) at
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
sed -i 's,${module_bedops},'"${module_bedops}"',g' "$logfile"
sed -i 's,${module_samtools},'"${module_samtools}"',g' "$logfile"
sed -i 's,${module_bedtools},'"${module_bedtools}"',g' "$logfile"
sed -i 's,${module_gmap},'"${module_gmap}"',g' "$logfile"
sed -i 's,${referencedir},'"${referencedir}"',g' "$logfile"
sed -i 's,${transcriptome},'"${transcriptome}"',g' "$logfile"
sed -i 's,${stem_genome},'"${stem_genome}"',g' "$logfile"
sed -i 's,${lncrnadb},'"${lncrnadb}"',g' "$logfile"
sed -i 's,${out_path_step16H1_BedTools},'"${out_path_step16H1_BedTools}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v
