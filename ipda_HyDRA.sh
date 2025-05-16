#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Mar 13, 2025
Last modified on May 16, 2025
Version: ${version}

#  _    _       _____  _____                    _            _ _            
# | |  | |     |  __ \|  __ \     /\           (_)          | (_)           
# | |__| |_   _| |  | | |__) |   /  \     _ __  _ _ __   ___| |_ _ __   ___ 
# |  __  | | | | |  | |  _  /   / /\ \   | '_ \| | '_ \ / _ \ | | '_ \ / _ \
# | |  | | |_| | |__| | | \ \  / ____ \  | |_) | | |_) |  __/ | | | | |  __/
# |_|  |_|\__, |_____/|_|  \_\/_/    \_\ | .__/|_| .__/ \___|_|_|_| |_|\___|
#          __/ |                         | |     | |                        
#         |___/                          |_|     |_|                        

Description: Execute individual HYDRA pipeline (HYbrid De novo RNA Assembly pipeline) command-lines for all steps.

Warning 1: PBS users can take advantage of the available write-to-pbs.sh scripts to write and submit PBS jobs for each step.
Warning 2: User must setup the settings (load/provide paths to all required tools) before executting each step.
Warning 3: Memory and CPU thresholds based on the HyDRA pipeline development are shown for each step, but these may vary across different datasets.
Warning 4: Do not run bash ipda_HYDRA.sh blindly:
            - Check log files and outputs at each step before proceeding.
            - Some steps (e.g., read processing) can run in parallel, e.g. short|long, but always verify step requirements before continuing.

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
"
}

#  _____ ______ _______    _____     _______ _    _  _____ 
# / ____|  ____|__   __|  |  __ \ /\|__   __| |  | |/ ____|
#| (___ | |__     | |     | |__) /  \  | |  | |__| | (___  
# \___ \|  __|    | |     |  ___/ /\ \ | |  |  __  |\___ \ 
# ____) | |____   | |     | |  / ____ \| |  | |  | |____) |
#|_____/|______|  |_|     |_| /_/    \_\_|  |_|  |_|_____/ 

# Set path to raw data and project
project_output_dir="/path/from/working/dir/to/Project_Name/"
short_raw_dir="/path/from/working/dir/to/data_short/" # Raw short_*1.f* and short_*2.f* files in fastq, fastq.gz, fq and fq.gz accepted
long_raw_dir="/path/from/working/dir/to/data_long/" # Raw long_read.f* file(s) in fastq, fastq.gz, fq and fq.gz accepted

# Set path to tools
# Disclaimer: Depending on how you have installed tools (See README.md file), you may have to provide path to other tools as well
fasta_splitter="/path/from/working/dir/to/fasta_splitter.py"
rcorrector="/path/from/working/dir/to/rcorrector/run_rcorrector.pl"
reformat="/path/from/working/dir/to/bbmap/reformat.sh"
fuper="/path/from/working/dir/to/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py"
bbduk="/path/from/working/dir/to/bbmap/bbduk.sh"
trinstats="/software/trinityrnaseq/trinityrnaseq-2.8.4/util/TrinityStats.pl"
picardtool="java -jar /software/picard/picard-tools-2.19.0/picard.jar"
gtf2gff3="/path/from/working/dir/to/gtf2gff3.pl"
genestats="docker/tools/genestats.sh"

# Set path to additional files
# WARNING: Description available next to each step
adapters_fasta="/path/from/working/dir/to/adapters_trimmomatic.fa"
ribosomal_rna_ref="/path/from/working/dir/to/ribosomalRNA.fa"
genome_bowtie2_index="/path/from/working/dir/to/ref_genome_bowtie2"
reference_gene="/path/from/working/dir/to/ref_gene.bed"
busco_odb9="/path/from/working/dir/to/BuscoLineages/eukaryota_odb9"
reference_fasta="/path/from/working/dir/to/ref-transcriptome.fasta"
reference_gmap_db=reference_GMAP_genome_db_dir
referencedir=/path/from/working/dir/to/
transcriptome_bed="/path/from/working/dir/to/ref-transcriptome.bed"
lncrnadb="/path/from/working/dir/to/lncRNAdb.bed"
lncrnadb_fasta="/path/from/working/dir/to/lncRNAdb.fasta" #bedtools getfasta -name -fi genome.fa -bed ${lncrnadb} > lncRNAdb.fasta
genome_proteins="/path/from/working/dir/to/protein-coding_reference-annotation.bed"
gencodegenome="/path/from/working/dir/to/gencode.genome.fa" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz>
gencodebed="/path/from/working/dir/to/gencode.v36.annotation.bed"
transgtf="/path/from/working/dir/to/ref-transcriptome.gtf"
ptngtf="/path/from/working/dir/to/gencode.v36.proteincoding.gtf" # Protein-coding transcripts were retrieved using grep "^#\|protein_coding"
lncrnagtf="/path/from/working/dir/to/gencode.v36.long_noncoding_RNAs.gtf" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.long_noncoding_RNAs.gtf.gz>
goodlncrnagtf="/path/from/working/dir/to/gencode.v36.confirmed_long_noncoding_RNAs.gtf" # grep -v "TEC" out of above (to be experimentally confirmed transcripts)


#  _____  ______          _____     _____  _____  ______      _____  _____   ____   _____ ______  _____ _____ _____ _   _  _____ 
# |  __ \|  ____|   /\   |  __ \|  |  __ \|  __ \|  ____|    |  __ \|  __ \ / __ \ / ____|  ____|/ ____/ ____|_   _| \ | |/ ____|
# | |__) | |__     /  \  | |  | |  | |__) | |__) | |__ ______| |__) | |__) | |  | | |    | |__  | (___| (___   | | |  \| | |  __ 
# |  _  /|  __|   / /\ \ | |  | |  |  ___/|  _  /|  __|______|  ___/|  _  /| |  | | |    |  __|  \___ \\___ \  | | | . ` | | |_ |
# | | \ \| |____ / ____ \| |__| |  | |    | | \ \| |____     | |    | | \ \| |__| | |____| |____ ____) |___) |_| |_| |\  | |__| |
# |_|  \_\______/_/    \_\_____/   |_|    |_|  \_\______|    |_|    |_|  \_\\____/ \_____|______|_____/_____/|_____|_| \_|\_____|


###############################################################
#                                                             #
#             STEP 01 [SHORT-READS]: QC RAW FILES             #
#                                                             #
###############################################################

## STEP 01S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
short_raw_1=${short_raw_dir}/short_read_1.fastq 
short_raw_2=${short_raw_dir}/short_read_2.fastq 
fastqc -t 1 --outdir ${project_output_dir}/hydra01S1_qc-raw_FastQC_${thislogdate} ${short_raw_1} ${short_raw_2} --memory 10000

## STEP 01S2 [resources:  -m 1 -c 1 -w "01:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
multiqc ${project_output_dir}/hydra01S1_qc-raw_FastQC_*/ -o ${project_output_dir}/hydra01S2_qc-raw_MultiQC_${thislogdate}

###############################################################
#                                                             #
#              STEP 01 [LONG-READS]: QC RAW FILES             #
#                                                             #
###############################################################

## STEP 01L1 [resources: -m 5 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw=${long_raw_dir}/long-read-stem.fastq 
long_raw_stem=long-read-stem
NanoPlot -t 1 -o ${project_output_dir}/hydra01L1_qc-raw_NanoPlot_${thislogdate}/${long_raw_stem} --fastq ${long_raw}

## STEP 01L2 [resources: -m 5 -c 1 -w "04:00:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw=${long_raw_dir}/long-read-stem.f* 
names="stem1 stem2"
NanoComp -t 1 --fastq ${long_raw} --names ${names} -o ${project_output_dir}/hydra01L2_qc-raw_NanoComp_${thislogdate}

## STEP 01L3 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw=${long_raw_dir}/long-read-stem.fastq 
fastqc -t 1 --outdir ${project_output_dir}/hydra01L3_qc-raw_FastQC_${thislogdate} ${long_raw} --memory 10000

## STEP 01L4 [resources: -m 1 -c 1 -w "00:30:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
multiqc ${project_output_dir}/hydra01L3_qc-raw_FastQC_* -o ${project_output_dir}/hydra01L4_qc-raw_MultiQC_${thislogdate}

## STEP 01L5 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw=${long_raw_dir}/long-read-stem.f*
ls ${long_raw} | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; seqtk seq -a ${file} > ${project_output_dir}/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta; done
ls ${long_raw} | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 ${fasta_splitter} -f ${project_output_dir}/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta  -n 20000 -p ${project_output_dir}/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}; done

###############################################################
#                                                             #
#        STEP 02 [SHORT-READS]: SHORT READS CORRECTION        #
#                                                             #
###############################################################

## STEP 02S1 [resources: -m 50 -c 15 -w "35:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
# WARNING 3: Set flag_merge to yes if you have samples sequenced from multiple lanes
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
short_raw_1=${short_raw_dir}/short_read_1.fastq 
short_raw_2=${short_raw_dir}/short_read_2.fastq 
flag_merge="yes"
if [[ $flag_merge = "yes" ]]; then zcat -c ${short_raw_1} | gzip -c >> ${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R1.fastq.gz ; fi
if [[ $flag_merge = "yes" ]]; then zcat -c ${short_raw_2} | gzip -c >> ${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R2.fastq.gz ; fi
if [[ $flag_merge = "yes" ]]; then f1="${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R1.fastq.gz"; f2="${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R2.fastq.gz"; else f1=${short_raw_1}; f2=${short_raw_2}; fi; perl ${rcorrector} -t 15 -od ${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate} -1 ${f1} -2 ${f2}
bash ${reformat} in=${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R1.cor.fq.gz in2=${project_output_dir}/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R2.cor.fq.gz

## STEP 02S2 [resources: -m 1 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
name="short_read_stem"
python ${fuper} -1 ${project_output_dir}/hydra02S1_short-cor_Rcorrector_DATE/short_read_R1.cor.fq.gz -2 ${project_output_dir}/hydra02S1_short-cor_Rcorrector_DATE/short_read_R2.cor.fq.gz -s ${name}
mv ${project_output_dir}/rmunfixable_${name}.log ${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/
mv ${project_output_dir}/unfixrm_${name}*1.cor.fq ${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/
mv ${project_output_dir}/unfixrm_${name}*2.cor.fq ${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/

###############################################################
#                                                             #
#          STEP 02 [LONG-READS]: LONG READS TRIMMING          #
#                                                             #
###############################################################

## STEP 02L1 [resources: -m 160 -c 40 -w "24:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw=${long_raw_dir}/long-read-stem.fastq 
long_raw_stem=long-read-stem
porechop -i ${long_raw} --threads 40 -o ${project_output_dir}/hydra02L1_trim-adapters_Porechop_${thislogdate}/adapter-trimmed_${long_raw_stem}.fastq.gz

## STEP 02L2 [resources: -m 2 -c 2 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw_stem=long-read-stem
gunzip -c ${project_output_dir}/hydra02L1_trim-adapters_Porechop_DATE/adapter-trimmed_${long_raw_stem}.fastq.gz | chopper -q 10 --headcrop 12 --threads 2 | gzip > ${project_output_dir}/hydra02L2_trim-qc-head_Chopper_${thislogdate}/quality-trimmed_adapter-trimmed_${long_raw_stem}.fastq.gz

## STEP 02L3 [resources: -m 1 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw_stem=long-read-stem
cutadapt --poly-a ${project_output_dir}/hydra02L2_trim-qc-head_Chopper_DATE/quality-trimmed_adapter-trimmed_${long_raw_stem}.fastq.gz -o ${project_output_dir}/hydra02L3_trim-polyA_Cutadapt_${thislogdate}/polyA-trimmed_quality-trimmed_adapter-trimmed_${long_raw_stem}.fastq.gz

###############################################################
#                                                             #
#          STEP 03 [SHORT-READS]: QC CORRECTED FILES          #
#                                                             #
###############################################################

## STEP 03S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
short_cor_1=${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_1.cor.fq 
short_cor_2=${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_2.cor.fq
fastqc -t 1 --outdir ${project_output_dir}/hydra03S1_qc-cor_FastQC_${thislogdate} ${short_cor_1} ${short_cor_2} --memory 10000

## STEP 03S2 [resources:  -m 1 -c 1 -w "01:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
multiqc ${project_output_dir}/hydra03S1_qc-cor_FastQC_*/ -o ${project_output_dir}/hydra03S2_qc-cor_MultiQC_${thislogdate}

###############################################################
#                                                             #
#           STEP 03 [LONG-READS]: QC TRIMMED FILES            #
#                                                             #
###############################################################

## STEP 03L1 [resources: -m 5 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_trim_stem=long_trim_stem
NanoPlot -t 1 -o ${project_output_dir}/hydra03L1_qc-trim_NanoPlot_${thislogdate}/${long_trim_stem} --fastq ${project_output_dir}/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fastq.gz

## STEP 03L2 [resources: -m 5 -c 1 -w "10:00:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
names="stem1 stem2"
NanoComp -t 1 --fastq ${project_output_dir}/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_*.fastq.gz --names ${names} -o ${project_output_dir}/hydra03L2_qc-trim_NanoComp_${thislogdate}

## STEP 03L3 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
fastqc -t 1 --outdir ${project_output_dir}/hydra03L3_qc-trim_FastQC_${thislogdate} ${project_output_dir}/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fastq.gz --memory 10000

## STEP 03L4 [resources: -m 1 -c 1 -w "00:30:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
multiqc ${project_output_dir}/hydra03L3_qc-trim_FastQC_DATE -o ${project_output_dir}/hydra03L4_qc-trim_MultiQC_${thislogdate}

## STEP 03L5 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
ls ${project_output_dir}/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_*.fastq.gz | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; seqtk seq -a ${file} > ${project_output_dir}/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta; done
ls ${long_raw} | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 ${fasta_splitter} -f ${project_output_dir}/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta  -n 20000 -p ${project_output_dir}/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}; done

###############################################################
#                                                             #
#         STEP 04 [SHORT-READS]: SHORT READS TRIMMING         #
#                                                             #
###############################################################

## STEP 04S1 [resources: -m 1 -c 6 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
# WARNING 2: User must provide path to adapter file. Read carefully:
# In order to define which adapter to trim (especially when sequencing
# details are not fully known) one can consult FastQC results and check
# for adapter contamination.
# These are examples of common adapters from 
# <http://docs.blast2go.com/user-manual/tools-(pro-feature)/fastq-quality-check/#FASTQQualityCheck-AdapterContent>
# Illumina Universal Adapter: AGATCGGAAGAG (12bp)
# Nextera Transposase Sequence: CTGTCTCTTATA (12bp)
adapters_fasta="/path/from/working/dir/to/adapters_trimmomatic.fa"
trimmomatic PE -phred33 -threads 6 ${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_stem_R1.cor.fq ${project_output_dir}/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_stem_R2.cor.fq ${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/adapter-trimmed_unfixrm_stem_R1.fq.gz ${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/unpaired_adapter-trimmed_unfixrm_stem_R1.fq.gz ${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/adapter-trimmed_unfixrm_stem_R2.fq.gz ${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/unpaired_adapter-trimmed_unfixrm_stem_R2.fq.gz ILLUMINACLIP:${adapters_fasta}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:75 HEADCROP:12

## STEP 04S2 [resources: -m 5 -c 15 -w "00:30:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
${bbduk} t=15 -Xmx5g in=${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R1.fq.gz in2=${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R2.fq.gz qtrim=rl trimq=30 minavgquality=20 minbasequality=10 minlen=75 outm=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-failed_stem_R1.cor.fq.gz outm2=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-failed_stem_R2.cor.fq.gz bhist=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/stem_bhist.txt qhist=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_qhist.txt gchist=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_gchist.txt aqhist=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_aqhist.txt lhist=${project_output_dir}/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_lhist.txt gcbins=auto

###############################################################
#                                                             #
#         STEP 04 [LONG-READS]: LONG READS CORRECTION         #
#                                                             #
###############################################################

## STEP 04L1 [resources: -m 60 -c 5 -w "10:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
gunzip -c ${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_* | awk 'NR % 4 == 2' | tr NT TN | ropebwt2 -LR | tr NT TN >> ${project_output_dir}/hydra04L1_short-reads-bwt_RopeBWT2_${thislogdate}/raw_all-short-reads.bwt

## STEP 04L2 [resources: -m 20 -c 2 -w "25:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
fmlrc2-convert -i ${project_output_dir}/hydra04L1_short-reads-bwt_RopeBWT2_DATE/raw_all-short-reads.bwt ${project_output_dir}/hydra04L2_short-reads-bwt_FMLRC2-convert_${thislogdate}/comp_msbwt_all-short-reads.npy

## STEP 04L3 [resources: -m 300 -c 8 -w "12:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
long_raw_stem=long-read-stem
fmlrc2 --threads 8 --cache_size 8 ${project_output_dir}/hydra04L2_short-reads-bwt_FMLRC2-convert_DATE/comp_msbwt_all-short-reads.npy ${project_output_dir}/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_${long_raw_stem}.fastq.gz ${project_output_dir}/hydra04L3_long-cor_FMLRC2_${thislogdate}/corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_${long_raw_stem}.fa
gzip ${project_output_dir}/hydra04L3_long-cor_FMLRC2_${thislogdate}/corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_${long_raw_stem}.fa

###############################################################
#                                                             #
#           STEP 05 [SHORT-READS]: QC TRIMMED FILES           #
#                                                             #
###############################################################

## STEP 05S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
short_trimm_1=${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_1.fq.gz 
short_trimm_2=${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_2.fq.gz
fastqc -t 1 --outdir ${project_output_dir}/hydra05S1_qc-trim_FastQC_${thislogdate} ${short_trimm_1} ${short_trimm_2} --memory 10000

## STEP 05S2 [resources:  -m 1 -c 1 -w "01:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
multiqc ${project_output_dir}/hydra05S1_qc-trim_FastQC_*/ -o ${project_output_dir}/hydra05S2_qc-trim_MultiQC_${thislogdate}

###############################################################
#                                                             #
#           STEP 05 [LONG-READS]: QC CORRECTED FILES          #
#                                                             #
###############################################################

## STEP 05L1 [resources: -m 1 -c 1 -w "00:30:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
ls ${project_output_dir}/hydra04L3_long-cor_FMLRC2_DATE/ | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; gunzip -c ${file}.fa.gz > ${project_output_dir}/hydra05L1_size-split_FastaSplitter_${thislogdate}/${main}.fasta; done
ls ${project_output_dir}/hydra05L1_size-split_FastaSplitter_${thislogdate}/${main}.fasta | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 ${fasta_splitter} -f ${file} -n 20000 -p ${project_output_dir}/hydra05L1_size-split_FastaSplitter_${thislogdate}/${main}

###############################################################
#                                                             #
#           STEP 06 [SHORT-READS]: FILTER rRNA READS          #
#                                                             #
###############################################################

## STEP 06S1 [resources: -m 2 -c 10 -w "08:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ribosomal RNAs file
# WARNING 3: Reason for fixing headers at the end of this step: Trimmomatic seems to change the naming scheme of reads. Original is "/1" for left read and "/2" for right reads respectively. Trinity needs the naming to be kept as expected. Source: https://www.biostars.org/p/141602/
# WARNING 4: Single End reads that survive trimming without a pair are discarded in the current pipeline version.
# WARNING 5: BBDuk uses more memory and CPUs than what is required. Therefore, we only ask for 60% of the memory and CPU specified when running this step.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
ribosomal_rna_ref="/path/from/working/dir/to/ribosomalRNA.fa"
${bbduk} t=6 -Xmx1g in=${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R1.fq.gz in2=${project_output_dir}/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R2.fq.gz outm=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribosomal-RNA_stem_R1.cor.fq.gz outm2=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribosomal-RNA_stem_R2.cor.fq.gz out=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R1.cor.fq out2=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R2.cor.fq ref=${ribosomal_rna_ref} stats=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/rRNA-alignment_stats.txt
cat ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R1.cor.fq | awk '{ if (NR%4==1) { print $1substr($2,2)"/1" } else { print } }' | sed 's/_1://g' >> ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R1.cor.fq; gzip ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R1.cor.fq; gzip ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R1.cor.fq
cat ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R2.cor.fq | awk '{ if (NR%4==1) { print $1substr($2,2)"/2" } else { print } }' | sed 's/_2://g' >> ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R2.cor.fq; gzip ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R2.cor.fq; gzip ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R2.cor.fq

###############################################################
#                                                             #
#           STEP 06 [LONG-READS]: FILTER rRNA READS           #
#                                                             #
###############################################################

## STEP 06L1 [resources: -m 5 -c 5 -w "00:30:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ribosomal RNAs file
# WARNING 3: BBDuk uses more memory and CPUs than what is required. Therefore, we only ask for 60% of the memory and CPU specified when running this step.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
ribosomal_rna_ref="/path/from/working/dir/to/ribosomalRNA.fa"
${bbduk} t=3 -Xmx3g in=${project_output_dir}/hydra04L3_long-cor_FMLRC2_DATE/corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fa.gz outm=${project_output_dir}/hydra06L1_ribodepletion_BBDuk_${thislogdate}/ribosomal-RNA_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fasta out=${project_output_dir}/hydra06L1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fasta ref=${ribosomal_rna_ref} stats=${project_output_dir}/hydra06L1_ribodepletion_BBDuk_${thislogdate}/rRNA-alignment_stats.txt

###############################################################
#                                                             #
#          STEP 07 [SHORT-READS]: MAP READS 2 GENOME          #
#                                                             #
###############################################################

## STEP 07S1 [resources: -m 30 -c 15 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to bowtie2 ref genome index
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
genome_bowtie2_index="/path/from/working/dir/to/ref_genome_bowtie2"
bowtie2 --fast --no-unal -p 15 -x ${genome_bowtie2_index} -1 ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R1.cor.fq.gz -2 ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R2.cor.fq.gz | samtools view -@15 -bS - >> ${project_output_dir}/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.bam
samtools sort -@15 ${project_output_dir}/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.bam -o ${project_output_dir}/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam
samtools index -@15 ${project_output_dir}/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam ${project_output_dir}/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam.bai

###############################################################
#                                                             #
#          STEP 07 [LONG-READS]: MAP READS 2 GENOME           #
#                                                             #
###############################################################

## STEP 07L1 [resources: -m 150 -c 15 -w "30:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to bowtie2 ref genome index
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
genome_bowtie2_index="/path/from/working/dir/to/ref_genome_bowtie2"
bowtie2 --fast --no-unal -p15 -x ${genome_bowtie2_index} -f ${project_output_dir}/hydra06L1_ribodepletion_BBDuk_DATE/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fasta | samtools view -@15 -bS - >> ${project_output_dir}/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.bam
samtools sort -@15 ${project_output_dir}/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.bam -o ${project_output_dir}/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.coordsorted.bam
samtools index -@15 ${project_output_dir}/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.bam ${project_output_dir}/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.coordsorted.bam.bai

###############################################################
#                                                             #
#          STEP 08 [SHORT-READS]: STRAND ASSESSMENT           #
#                                                             #
###############################################################

## STEP 08S1 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ref gene BED file. e.g. download <https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz>, unzip files and use sed to change chromosome nomenclature (to match your reference genome): sed 's/^chr//g'. WARNING! save files as hg38RefGene.bed
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
reference_gene="/path/from/working/dir/to/ref_gene.bed"
infer_experiment.py -r ${reference_gene} -i ${project_output_dir}/hydra07S1_gen-map_Bowtie2_DATE/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam >> ${project_output_dir}/hydra08S1_strandness-assess_RSeQC_${thislogdate}/strand-assessment_genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.txt

## STEP 08S2 [resources: -m 1 -c 8 -w "02:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
hex_name_colors="black;grey;lightgrey;palevioletred;firebrick;crimson;orangered;saddlebrown;bisque;navajowhite;orchid;magenta;darkkhaki;olive;darkolivegreen;darkseagreen;darkgreen;mediumseagreen;aquamarine;lightcyan;darkcyan;cadetblue;lightskyblue;lightslategrey;royalblue;darkblue;darkslateblue;indigo;plum;fuchsia;hotpink;lightpink" # if you have more than 30 samples you will need to add more colours to this list
#hex_name_colors="#61bcc6;#ade66c;#f11e71;#9cd5d8;#6213fe;#e17ebb;#f7dc02;#a14b74;#a66fc7;#de08e3;#be1f5f;#1838f9;#d6fd2e;#16a184;#5b1119;#b9f9f7;#57e148;#0c4b0e;#a88f99;#4b8623;#0c170d;#cfc154;#424cab;#669dd7;#3a06a8;#ec5a5e;#56dac1;#e4d2d6;#262f4f;#b3d17e"
colors=`echo ${hex_name_colors} | tr ';' '\n' | tr '\n' ' ' | head -n${sample_num}`
labels="stem1 stem2"
echo '${hex_name_colors}' | tr ';' '\n' >> ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-pallete_${thislogdate}.txt
echo '${labels}' | tr ' ' '\n' | uniq >> ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels_${thislogdate}.txt
head -n${sample_num} ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-pallete_${thislogdate}.txt >> ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-scheme_${thislogdate}.txt
paste ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels_${thislogdate}.txt ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-scheme_${thislogdate}.txt >> ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels-color-scheme_${thislogdate}.txt
cut -f2 ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels-color-scheme_${thislogdate}.txt >> ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-sample-colors_${thislogdate}.txt
multiBamSummary bins --bamfiles ${project_output_dir}/hydra07S1_gen-map_Bowtie2_DATE/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm*.coordsorted.bam --labels ${labels} --outFileName ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --maxFragmentLength 280 --numberOfProcessors 8
plotCorrelation --corData ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --corMethod spearman --whatToPlot heatmap --labels ${labels} --plotFile ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/spearman-correlation_analysis-stem.pdf --plotTitle Spearman-Correlation_analysis-stem_${thislogdate} --outFileCorMatrix ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/spearman-correlation_analysis-stem.matrix --plotNumbers
plotCorrelation --corData ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --corMethod pearson --whatToPlot heatmap --labels ${labels} --plotFile ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pearson-correlation_analysis-stem.pdf --plotTitle Pearson-Correlation_analysis-stem_${thislogdate} --outFileCorMatrix ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pearson-correlation_analysis-stem.matrix --plotNumbers
plotPCA --corData ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --plotHeight 20 --labels ${labels} --colors ${colors} --plotFile ${project_output_dir}/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-plot_analysis-stem.pdf --plotTitle PCA-plot_analysis-stem_${thislogdate} --log2

###############################################################
#                                                             #
#           STEP 08 [LONG-READS]: STRAND ASSESSMENT           #
#                                                             #
###############################################################

## STEP 08L1 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ref gene BED file. e.g. download <https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz>, unzip files and use sed to change chromosome nomenclature (to match your reference genome): sed 's/^chr//g'. WARNING! save files as hg38RefGene.bed
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
reference_gene="/path/from/working/dir/to/ref_gene.bed"
infer_experiment.py -r ${reference_gene} -i ${project_output_dir}/hydra07L1_gen-map_Bowtie2_DATE/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.coordsorted.bam >> ${project_output_dir}/hydra08L1_strandness-assess_${thislogdate}/strand-assessment_genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.txt

###############################################################
#                                                             #
#        STEP 09 [SHORT-READS]: QC RIBODEPLETED FILES         #
#                                                             #
###############################################################

## STEP 09S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
short_ribodepleted_1=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_1.cor.fq.gz 
short_ribodepleted_2=${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_2.cor.fq.gz 
fastqc -t 1 --outdir ${project_output_dir}/hydra09S1_qc-ribodepleted_FastQC_${thislogdate} ${short_ribodepleted_1} ${short_ribodepleted_2} --memory 10000

## STEP 09S2 [resources:  -m 1 -c 1 -w "01:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
multiqc ${project_output_dir}/hydra09S1_qc-ribodepleted_FastQC_*/ -o ${project_output_dir}/hydra09S2_qc-ribodepleted_MultiQC_${thislogdate}

###############################################################
#                                                             #
#         STEP 09 [LONG-READS]: QC RIBODEPLETED FILES         #
#                                                             #
###############################################################

## STEP 09L1 [resources: -m 1 -c 1 -w "00:30:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
ls ${project_output_dir}/hydra06L1_ribodepletion_BBDuk_DATE/${main}.fasta | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 ${fasta_splitter} -f ${file} -n 20000 -p ${project_output_dir}/hydra09L1_size-split_FastaSplitter_${thislogdate}/${main}

###############################################################
#                                                             #
#    STEP 10 [HYBRID]: SUMMARY QC AND FILTERING RAW READS     #
#                                                             #
###############################################################

## STEP 10H1 [resources: -m 1 -c 1 -w "24:00:00"]
# This step produces a series of TSV files from out and error files printed from stdout on a PBS. Please adapt where necessary to suit how you've run the steps so far.
# Please go through each  series of commands on write-to-pbs/ipda_HYDRA_step10H1-to-pbs.sh

# _______                            _       _                                                       _     _       
#|__   __|                          (_)     | |                             /\                          | |   | |      
#   | |_ __ __ _ _ __  ___  ___ _ __ _ _ __ | |_ ___  _ __ ___   ___       /  \   ___ ___  ___ _ __ ___ | |__ | |_   _ 
#   | | '__/ _` | '_ \/ __|/ __| '__| | '_ \| __/ _ \| '_ ` _ \ / _ \     / /\ \ / __/ __|/ _ \ '_ ` _ \| '_ \| | | | |
#   | | | | (_| | | | \__ \ (__| |  | | |_) | || (_) | | | | | |  __/    / ____ \\__ \__ \  __/ | | | | | |_) | | |_| |
#   |_|_|  \__,_|_| |_|___/\___|_|  |_| .__/ \__\___/|_| |_| |_|\___/  _/    \_\___/___/\___|_| |_| |_|_.__/|_|\__, |
#                                     | |                                                                     __/ |
#                                     |_|                                                                    |___/ 
#

###############################################################
#                                                             #
#      STEP 11 [HYBRID]: DE NOVO TRANSCRIPTOME ASSEMBLY       #
#                                                             #
###############################################################

## STEP 11H1 [resources: -m 250 -c 40 -w "60:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
zcat -c ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm*_R1.cor.fq.gz | gzip -c >> ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R1.cor.fq.gz
zcat -c ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm*_R2.cor.fq.gz | gzip -c ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R2.cor.fq.gz
gzip -c ${project_output_dir}/hydra06L1_ribodepletion_BBDuk_DATE/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_*.fasta >> ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-long-reads_merged.fa.gz
rnaspades.py --checkpoint "all" -t 40 -m 250 -1 ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R1.cor.fq.gz -2 ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R2.cor.fq.gz --nanopore ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-long-reads_merged.fa.gz -o ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}
seqtk seq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/transcripts.fasta > ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/transcripts.singleline.fasta
rm -f ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R1.cor.fq.gz ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R2.cor.fq.gz ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-long-reads_merged.fa.gz

## STEP 11H2 [resources: -m 250 -c 12 -w "85:00:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
ls ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R1.cor.fq.gz >> trinity-left-list_${thislogdate}.txt
ls ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R2.cor.fq.gz >> trinity-right-list_${thislogdate}.txt
left_list=`cat trinity-left-list_${thislogdate}.txt | tr '\n' ',' | tr -d ' '`; left_list="${left_list%\,}"; right_list=`cat trinity-right-list_${thislogdate}.txt | tr '\n' ',' | tr -d ' '`; right_list="${right_list%\,}"
Trinity --seqType fq --max_memory 250G --left $left_list --right $right_list --CPU 12 --SS_lib_type RF --output ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate} --group_pairs_distance 280 --full_cleanup --bflyHeapSpaceMax 16G --normalize_max_read_cov 50 --bypass_java_version_check
mv trinity-left-list_${thislogdate}.txt trinity-right-list_${thislogdate}.txt ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}
mv ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}.Trinity.fasta  ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}/assembly.Trinity.fasta
mv ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}.Trinity.fasta.gene_trans_map ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}/assembly.Trinity.fasta.gene_trans_map

###############################################################
#                                                             #
#          STEP 12 [HYBRID]: QUALITY CHECK ASSEMBLY           #
#                                                             #
###############################################################

## STEP 12H1 [resources: -m 2 -c 5 -w "50:00:00"]
# WARNING: User must provide path to BUSCO lineage dataset, which can be obtained from <https://busco-archive.ezlab.org/v3/datasets/eukaryota_odb9.tar.gz>. Filename eukaryota_odb9 given as an example
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
busco_odb9="/path/from/working/dir/to/BuscoLineages/eukaryota_odb9"
tar -xvf /software/busco/busco-20161219/augustus-3.2.2.config.tgz # only if not previously done, untar augustus in the home folder (outside Interactive session in a PBS HPC)
BUSCO.py -i ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta -l ${busco_odb9} -o eukaryota-odb9_stem-rnaSPAdes -m transcriptome -c 5 -t ${project_output_dir}/hydra12H1_qc_BUSCO_${thislogdate} -sp human
BUSCO.py -i ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta -l ${busco_odb9} -o eukaryota-odb9_stem-Trinity -m transcriptome -c 5 -t ${project_output_dir}/hydra12H1_qc_BUSCO_${thislogdate} -sp human

## STEP 12H2 [resources: -m 1 -c 1 -w "00:30:00"]
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
${trinstats} ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta >> ${project_output_dir}/hydra12H2_qc_TrinityStats_${thislogdate}/trinityStats_stem-rnaSPAdes
${trinstats} ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta >> ${project_output_dir}/hydra12H2_qc_TrinityStats_${thislogdate}/trinityStats_stem-Trinity

## STEP 12H3 [resources: -m 80 -c 15 -w "120:00:00"]
# WARNING: User must provide path to reference transcriptome annotation file
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
reference_fasta="/path/from/working/dir/to/ref-transcriptome.fasta"
transrate --assembly ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta --reference ${reference_fasta} --threads 15 --output ${project_output_dir}/hydra12H3_qc_TransRate_${thislogdate}/
mv ${project_output_dir}/hydra12H3_qc_TransRate_${thislogdate}/assemblies.csv ${project_output_dir}/hydra12H3_qc_TransRate_${thislogdate}/stem_assemblies.csv

## STEP 12H4 [resources: -m 85 -c 32 -w "05:00:00"]
# WARNING: User must provide path to reference transcriptome annotation file
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
reference_fasta="/path/from/working/dir/to/ref-transcriptome.fasta"
tablesize=`awk 'END{print NR}' hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq | cut -d' ' -f1`
fastq_pair -t ${tablesize} hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/right.norm.fq
transrate --assembly ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta  --left hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq --right hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/right.norm.fq.paired.fq --reference ${reference_fasta} --threads 32 --output ${project_output_dir}/hydra12H4_qc-and-polish_TransRate_${thislogdate}
mv ${project_output_dir}/hydra12H4_qc-and-polish_TransRate_${thislogdate}/assemblies.csv ${project_output_dir}/hydra12H4_qc-and-polish_TransRate_${thislogdate}/stem_assemblies.csv

###############################################################
#                                                             #
#     STEP 13 [HYBRID]: TRANSCRIPTOME READ REPRESENTATION     #
#                                                             #
###############################################################

## STEP 13H1 [resources: -m 250 -c 16 -w "120:00:00"]
# WARNING: You will need to loop through sample files for Bowtie2 and GMAP alignments.
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
bowtie2-build --threads 16 ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_bowtie2idx
bowtie2 -p 16 -q --no-unal -k 20 -x ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_bowtie2idx -1 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R1.cor.fq.gz -2 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R2.cor.fq.gz | samtools view -@16 -Sb -o ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.bam
gmap -n 0 -t 16 -B 5 --dir=${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/ --gseg ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta ${project_output_dir}/hydra06L1_ribodepletion_BBDuk_DATE/stem.fasta | samtools view -@16 -bS - >> ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.bam
samtools merge -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.bam ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.bam ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.bam
samtools sort -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.bam -o ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.coordsorted.bam
samtools sort -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.bam -o ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam
samtools sort -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.bam -o ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.coordsorted.bam
samtools index -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.coordsorted.bam
samtools index -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam
samtools index -@16 ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.coordsorted.bam
samtools faidx ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta -o ${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_samtoolsidx
${picardtool} CollectInsertSizeMetrics I=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.coordsorted.bam O=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-transcriptome-alignment-merged.txt H=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-merged-histogram.pdf INCLUDE_DUPLICATES=FALSE
${picardtool} CollectInsertSizeMetrics I=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam O=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-transcriptome-alignment-shortreads.txt H=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-shortreads-histogram.pdf INCLUDE_DUPLICATES=FALSE
${picardtool} CollectInsertSizeMetrics I=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.coordsorted.bam O=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-transcriptome-alignment-longreads.txt H=${project_output_dir}/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-longreads-histogram.pdf INCLUDE_DUPLICATES=FALSE

## STEP 13H2 [resources: -m 15 -c 15 -w "05:00:00"]
# WARNING: You will need to loop through sample files for Bowtie2 alignments.
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
bowtie2-build --threads 15 ${project_output_dir}/${project_output_dir}/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_bowtie2idx
bowtie2 -p 15 -q --no-unal -k 20 -x ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_bowtie2idx -1 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R1.cor.fq.gz -2 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R2.cor.fq.gz | samtools view -@15 -Sb -o ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.bam
samtools sort -@15 ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.bam -o ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam
samtools index -@15 ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam
samtools faidx ${project_output_dir}/${project_output_dir}/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta -o ${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_samtoolsidx
${picardtool} CollectInsertSizeMetrics I=${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-longreads-gmap.coordsorted.bam O=${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_insert-size-transcriptome-alignment-longreads.txt H=${project_output_dir}/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_insert-size-longreads-histogram.pdf INCLUDE_DUPLICATES=FALSE

###############################################################
#                                                             #
#            STEP 14 [HYBRID]: SPLICING ASSESSMENT            #
#                                                             #
###############################################################

## STEP 14H1 [resources: -m 15 -c 16 -w "10:00:00"]
# WARNING 1: User must provide name of reference_GMAP_genome_db_dir and its path separately.
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
reference_gmap_db=reference_GMAP_genome_db_dir
referencedir=/path/from/working/dir/to/
gmap -n 0 -t 16 -B 5 --dir=${referencedir} --db=${reference_gmap_db} ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts_splicing-assessment.genome.sam
grep -v "^@" ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts_splicing-assessment.genome.sam | awk -v A="${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.duplicated-alignments.sam" -v B="${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam" '{count[$1]++; lines[$1,count[$1]] = $0;} END {for (key in count) {file = (count[key] == 2 ? A : B); for (i = 1; i <= count[key]; i++) print lines[key, i] > file;}}'
sed -n 'p;n' ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.duplicated-alignments.sam >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.sam
cat ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam | awk '$1 ~ /@/ || ($6 !~ /N/ && $6 !~ /\*/)' >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.sam
cat ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam | awk '$1 ~ /@/ || ($6 ~ /N/)' >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.sam
cat ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam | awk '$1 ~ /@/ || ($6 ~ /\*/)' >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unaligned.sam
grep -v "^@" ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.sam | cut -f1 >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.ids
grep -v "^@" ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.sam | cut -f1 >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.ids
seqtk subseq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.ids >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.fasta
seqtk subseq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.ids >> ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.fasta

###############################################################
#                                                             #
#                STEP 15 [HYBRID]: READ SUPORT                #
#                                                             #
###############################################################

## STEP 15H1 [resources: -m 20 -c 16 -w "20:00:00"]
# WARNING: Reproduce for unspliced
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
identity=98; coverage=90 # CDhit cutoffs
high_rpkm_cutoff=1 # e.g.: 1 for spliced and 5 for unspliced
low_rpkm_cutoff=0.5 # e.g.: 0.5 for spliced and 3 for unspliced
cd-hit-est -o ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta -c 0.${identity} -i ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_DATE/transcripts/transcripts.spliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T 16 -M 20000
minimap2 -ax splice -t 16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra06L1_ribodepletion_BBDuk_DATE/stem1.fasta,${project_output_dir}/hydra06L1_ribodepletion_BBDuk_DATE/stem2.fasta | samtools view -@16 -bS - >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.bam
bowtie2-build --threads 16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx
bowtie2 -q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --threads 16 -x ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx -1 ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq -2 ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq | samtools view -@16 -bS - >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam
samtools merge ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap.spliced.bam -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.bam
samtools sort -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam -o ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.coordsorted.bam
samtools sort -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.bam -o ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.coordsorted.bam
samtools sort -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap.spliced.bam -o ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap.spliced.coordsorted.bam
rm -f ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.bam ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap.spliced.bam
samtools index -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.coordsorted.bam
samtools index -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.coordsorted.bam
samtools index -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap.spliced.coordsorted.bam
samtools idxstats -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.coordsorted.bam >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts.spliced.tsv
samtools idxstats -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap.spliced.coordsorted.bam >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap_alignment-counts.spliced.tsv
samtools idxstats -@16 ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap.spliced.coordsorted.bam >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap_alignment-counts.spliced.tsv
total=`sed '$d' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv
total=`sed '$d' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv
total=`sed '$d' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv
awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta
awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta
awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta
awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta
awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_longreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta
awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_shortreads_bt2_and_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_merged_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta

## STEP 15H2 [resources: -m 20 -c 16 -w "20:00:00"]
# WARNING: Reproduce for unspliced
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
identity=98; coverage=90 # CDhit cutoffs
high_rpkm_cutoff=1 # e.g.: 1 for spliced and 5 for unspliced
low_rpkm_cutoff=0.5 # e.g.: 0.5 for spliced and 3 for unspliced
cd-hit-est -o ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta -c 0.${identity} -i ${project_output_dir}/hydra14H1_splicing-assessment_GMAP_DATE/transcripts/transcripts.spliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.${coverage} -p 1 -d 0 -b 3 -T 16 -M 20000
bowtie2-build --threads 16 ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx
bowtie2 -q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --threads 16 -x ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx -1 ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq -2 ${project_output_dir}/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq | samtools view -@16 -bS - >> ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam
samtools sort -@16 ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam -o ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.coordsorted.bam
rm -f ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.bam
samtools index -@16 ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.coordsorted.bam
samtools idxstats -@16 ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2.spliced.coordsorted.bam >> ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts.spliced.tsv
total=`sed '$d' ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv
awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta
awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}.spliced.fasta ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> ${project_output_dir}/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity${identity}_cov${coverage}_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta

###############################################################
#                                                             #
#         STEP 16 [HYBRID]: TRANSCRIPTOME ANNOTATION          #
#                                                             #
###############################################################

## STEP 16H1 [resources: -m 20 -c 16 -w "10:00:00"]
# WARNING 1: User must provide: name of reference_GMAP_genome_db_dir and its path separately; reference transcriptome BED file; lncRNA database BED file
# WARNING 2: Reproduce for unspliced
# WARNING 3: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
reference_gmap_db=reference_GMAP_genome_db_dir
referencedir=/path/from/working/dir/to/
stem_genome="refgenome_stem"
transcriptome_bed="/path/from/working/dir/to/ref-transcriptome.bed"
lncrnadb="/path/from/working/dir/to/lncRNAdb.bed"
gmap -n 1 -t 16 -B 5 --dir=${referencedir} --db=${reference_gmap_db} ${project_output_dir}/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse > ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.sam
convert2bed-megarow --input=sam --do-not-sort --max-mem=20G < ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.sam > ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.bed12
cut -f1-6 ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.bed12 >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.bed
bedtools intersect -wo -s -f 0.75 -r -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.bed -b ${transcriptome} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomefeatures.bed
bedtools intersect -s -f 0.75 -r -C -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.bed -b ${transcriptome} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomecounts.bed
bedtools intersect -v -wa -s -f 0.75 -r -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.bed -b ${transcriptome} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.notingenome.bed
bedtools intersect -wo -s -f 0.75 -r -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.notingenome.bed -b ${lncrnadb} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.lncRNAfeatures.bed
bedtools intersect -s -f 0.75 -r -C -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.notingenome.bed -b ${lncrnadb} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.lncRNAcounts.bed
bedtools intersect -v -wa -s -f 0.75 -r -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.notingenome.bed -b ${lncrnadb} ${transcriptome} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.novel.bed
bedtools intersect  -wo -s -f 0.5 -e -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.novel.bed -b ${lncrnadb} ${transcriptome} >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.potentialfragment.bed
cut -f4 ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomefeatures.bed >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_a ; cut -d'"' -f2,6 ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomefeatures.bed | tr '"' '\t' >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_b ; grep -Po 'gene_biotype "\K.*?(?=")' ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomefeatures.bed >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_c ; grep -Po 'transcript_biotype "\K.*?(?=")' ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomefeatures.bed >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_d ; paste ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_b ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_c ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_d >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation; rm -rf ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_b ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_c ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.tmp_d ; cut -f4,10,11 ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.lncRNAfeatures.bed | sed "s/\t./\.\tncRNAdb\t\./g" >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation
cut -f1,4 ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation | sort | uniq >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary
grep "protein_coding" ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary.protein.transcriptID
grep "lincRNA\|ncRNAdb\|lncRNA" ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary.ncRNA.transcriptID
grep "antisense" ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary.antisense.transcriptID
grep "pseudogene" ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.annotation.summary.pseudogene.transcriptID
cut -f12,16 ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}.genomefeatures.bed | tr '"' '\t' | sort | uniq >> ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.${stem_genome}_geneID-type_retrieved-annotated-genes

###############################################################
#                                                             #
#    STEP 17 [HYBRID]: SUMMARY ASEEMBLY AND QC RAW READS     #
#                                                             #
###############################################################

## STEP 17H1 [resources: -m 1 -c 1 -w "24:00:00"]
# This step produces a series of TSV files from out and error files printed from stdout on a PBS. Please adapt where necessary to suit how you've run the steps so far.
# Please go through each  series of commands on write-to-pbs/ipda_HYDRA_step17H1-to-pbs.sh

# _            _____  _   _                _ #_                                   
#| |          |  __ \| \ | |   /\         | (_)                                  
#| |_ __   ___| |__) |  \| |  /  \      __| |_ ___  ___ _____   _____ _ __ _   _ 
#| | '_ \ / __|  _  /| . ` | / /\ \    / _` | / __|/ __/ _ \ \ / / _ \ '__| | | |
#| | | | | (__| | \ \| |\  |/ ____ \   | (_| | \__ \ (_| (_) \ V /  __/ |  | |_| |
#|_|_| |_|\___|_|  \_\_| \_/_/    \_\  \__,_|_|___/\___\___/ \_/ \___|_|   \__, |
#                                                                        __/ |
#                                                                       |___/ 

###############################################################
#                                                             #
#         STEP 18 [HYBRID]: PREDICT CODING POTENTIAL          #
#                                                             #
###############################################################

## STEP 18H1 [resources: -m 10 -c 1 -w "400:00:00"]
# WARNING 1: PLEK may have switched the labels "coding" and "non-coding"
# WARNING 2: Reproduce for unspliced
# WARNING 3: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
genome_proteins="/path/from/working/dir/to/protein-coding_reference-annotation.bed"
bedtools intersect -wo -s -f 0.6 -r -a ${project_output_dir}/hydra16H1_transcriptome-annotation_BedTools_DATE/transcripts/transcripts.bed -b ${genome_proteins} >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.proteincoding.bed
ezLncPred -i ${project_output_dir}/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -o ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2 CPC2
ezLncPred -i ${project_output_dir}/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -o ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT CPAT -p Human
ezLncPred -i ${project_output_dir}/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -o ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI CNCI -p ve
cut -f1,7,8 ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2 | grep -w "coding" | sort | uniq >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.proteincoding
cut -f1,7,8 ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2 | grep -w "noncoding" | sort | uniq >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.noncoding
cat ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT | tr "'" '\t' | cut -f2,8 | awk 'NR>0{if($2>=0.5){print}}' | sort | uniq >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.proteincoding
cat ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT | tr "'" '\t' | cut -f2,8 | awk 'NR>0{if($2<0.5){print}}' | sort | uniq >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.noncoding
cat ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI/CNCI.index | tail -n +2 | cut -f1,2,3 | grep -Pe "\tcoding" | awk 'NR>0{printf "%s\t%.2f\n", $1,$3}' | sort | uniq >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI.proteincoding
cat ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI/CNCI.index | tail -n +2 | cut -f1,2,3 | grep -Pe "\tnoncoding" | awk 'NR>0{printf "%s\t%.2f\n", $1,$3}' | sort | uniq >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI.noncoding
number=`cut -f1 ${annotation} ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.proteincoding | tail -n +2 | sort | uniq -c | grep -c " 2 "` ; echo "${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts\tCPC2 Annotated ptcoding classified as ptcoding\t${number}" >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.proteincodingtest
number=`cut -f1 ${annotation} ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.proteincoding | tail -n +2 | sort | uniq -c | grep -c " 2 "` ; echo "${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts\tCPAT Annotated ptcoding classified as ptcoding\t${number}" >> ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.proteincodingtest

###############################################################
#                                                             #
#         STEP 19 [HYBRID]: IDENTIFY lncRNA POTENTIAL         #
#                                                             #
###############################################################

## STEP 19H1 [resources: -m 120 -c 16 -w "60:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
gencodegenome="/path/from/working/dir/to/gencode.genome.fa" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz>
gencodebed="/path/from/working/dir/to/gencode.v36.annotation.bed"
transgtf="/path/from/working/dir/to/ref-transcriptome.gtf"
reference_fasta="/path/from/working/dir/to/ref-transcriptome.fasta"
ptngtf="/path/from/working/dir/to/gencode.v36.proteincoding.gtf" # Protein-coding transcripts were retrieved using grep "^#\|protein_coding"
lncrnagtf="/path/from/working/dir/to/gencode.v36.long_noncoding_RNAs.gtf" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.long_noncoding_RNAs.gtf.gz>
goodlncrnagtf="/path/from/working/dir/to/gencode.v36.confirmed_long_noncoding_RNAs.gtf" # grep -v "TEC" out of above (to be experimentally confirmed transcripts)
minimap2 -x splice --secondary=no -a -t 16 -o ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sam ${gencodegenome} ${project_output_dir}/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta
samtools view -u ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sam | samtools sort -@ 16 -o ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sorted.bam
bamToBed -splitD -bed12 -i ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sorted.bam > ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.bed12
bedtools intersect -wo -s -f 0.75 -r -a ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.bed12 -b ${gencodebed} >> ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.gencode.genomefeatures.bed
bedToGenePred ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.bed12 stdout | genePredToGtf file stdin ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.gtf
FEELnc_filter.pl -i ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.gtf -a $transgtf --biotype transcript_biotype=protein_coding --size=200 --monoex=0 -p 16 -o ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.feelncfilter.log > ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.lncRNAcandidates.gtf
FEELnc_codpot.pl -i ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.lncRNAcandidates.gtf -a $ptngtf -l $goodlncrnagtf -g $gencodegenome -o transcripts --outdir=${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts
FEELnc_classifier.pl -i ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/transcripts.lncRNA.gtf -a $ptngtf -l ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/transcripts.classifier.log >> ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/transcripts.lncRNAclasses.txt
mkdir -p ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/VarImpPlots ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/ROCcurves
mv ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/*varImpPlot.png ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/VarImpPlots ; mv ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/*TGROC.png ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/ROCcurves

###############################################################
#                                                             #
#              STEP 20 [HYBRID]: DEFINE lncRNAs               #
#                                                             #
###############################################################

## STEP 20H1 [resources: -m 1 -c 1 -w "01:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
grep -P "^1\t" ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt | cut -f3- >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.noncoding
cut -f1 ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_DATE/stem/stem.CNCI.noncoding ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_DATE/stem/stem.CPAT.noncoding ${project_output_dir}/hydra18H1_coding-potential_ezLncPred_DATE/stem/stem.CPC2.noncoding | sort | uniq -i -c | grep -v " 1 " | cut -c9- >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred.min2-models.noncoding
cut -f1 ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred.min2-models.noncoding ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.noncoding | sort | uniq -i -c | grep " 2 " | cut -c9- >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred-and-feelnc.noncoding
grep "isBest\|^1" ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.best; head -n1 ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered; head -n1 ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best; head -n1 ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best.intergenic; awk 'NR==FNR{a[tolower($1)]; next} tolower($2) in a' ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred-and-feelnc.noncoding ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.best >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best; cat ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best | awk 'NR==1{print}; NR>1{if($7=="intergenic"){print}}' >> ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best.intergenic

###############################################################
#                                                             #
#          STEP 21 [HYBRID]: RETRIEVE lncRNA METRICS          #
#                                                             #
###############################################################

## STEP 21H1 [resources: -m 1 -c 3 -w "01:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
cat ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_DATE/transcripts_spliced/transcripts_spliced.ezlncpred-and-feelnc.noncoding | while read transcript; do grep "gene_id \"${transcript}\";" ${project_output_dir}/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts_spliced/transcripts_spliced.gtf >> ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gtf; done
perl ${gtf2gff3} ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gtf >> ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gff
${genestats} ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gff >> ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.genestats
cat ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.genestats | while read line; do echo "${line}" | awk -v OFS="\t" '{ exon_sum += $3; exon_len += $4; intron_sum += $5; intron_len += $6 } END { print exon_sum, exon_len / exon_sum, intron_sum, intron_len / (intron_sum+1) }' >> ${project_output_dir}/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.metrics ; done

## STEP 21H2 [resources: -m 2 -c 20 -w "60:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
transcriptome_bed="/path/from/working/dir/to/ref-transcriptome.bed"
reference_fasta="/path/from/working/dir/to/ref-transcriptome.fasta"
transgtf="/path/from/working/dir/to/ref-transcriptome.gtf"
lncrnadb="/path/from/working/dir/to/lncRNAdb.bed"
lncrnadb_fasta="/path/from/working/dir/to/lncRNAdb.fasta" #bedtools getfasta -name -fi genome.fa -bed ${lncrnadb} > lncRNAdb.fasta
cut -f2 ${project_output_dir}/hydra20H1_predicted-lncRNAs_Bash_DATE/transcripts_spliced/transcripts.spliced.feelncmain.filtered-and-best | sort | uniq >> ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.ids ; seqtk subseq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.singleline.fasta ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.ids > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.fasta
pblat -threads=20 -minIdentity=75 ${reference_fasta} ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.fasta ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.psl
pblat -threads=20 -minIdentity=75 ${lncrnadb_fasta} ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.fasta ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.psl
awk -F'\t' '{if($1>=0.85*$11 && $1>=0.85*$15){print}}' ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.psl | sed -e '1,2d' > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.psl
cut -f14 ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.psl | while read t; do grep -P "\ttranscript\t.*${t}" ${transgtf} ; done | sort | uniq > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.gtf
cat ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.psl | while read l; do t=`echo ${l} | cut -d" " -f14`; type=`grep "${t}" ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.gtf | cut -d'"' -f20` ; echo -e "${l}\t${type}"; done > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.types
awk -F'\t' '{if($1>=0.85*$11 && $1>=0.85*$15){print}}' ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.psl | sed -e '1,2d' > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.85c.psl
cat ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.85c.psl | while read l; do echo -e "${l}\tlncrnadb"; done >> ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.types
cat ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.types | while read l; do type=`echo ${l} | cut -d" " -f22`; id=`echo ${l} | cut -d" " -f10`; if [[ "$type" == "lincRNA" || "$type" == "lncrnadb" ]]; then echo -e "${id}" >> ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.lnc; else echo -e "${id}" >> ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc.tmp; fi; done ; sort ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc.tmp | uniq > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc; rm -f ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc.tmp
sort ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.ids ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc | uniq -c | grep " 1 " | tr -s ' ' | cut -d' ' -f3 > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc
cat ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc | while read t; do grep -P "\ttranscript\t.*${t}" ${transgtf} ; done > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc.gtf
seqtk subseq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.singleline.fasta ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc > ${project_output_dir}/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc.fasta

###############################################################
#                                                             #
#         STEP 22 [HYBRID]: SUMMARY lncRNA DISCOVERY          #
#                                                             #
###############################################################

## STEP 22H1 [resources: -m 1 -c 1 -w "24:00:00"]
# This step produces a series of TSV files from out and error files printed from stdout on a PBS. Please adapt where necessary to suit how you've run the steps so far.
# Please go through each  series of commands on write-to-pbs/ipda_HYDRA_step22H1-to-pbs.sh 
