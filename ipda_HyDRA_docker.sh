#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Apr 24, 2025
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

Description: Execute individual HYDRA pipeline (HYbrid De novo RNA Assembly pipeline) command-lines for all steps using the Docker container.

Warning 1: HyDRA is also available in step by step BASH command-lines and write-to-pbs.sh BASH scripts. NextFlow integration coming!
Warning 2: User must updated /path/from/root/to/project/data/ before executting each step.
Warning 3: Memory and CPU thresholds are based on the HyDRA pipeline development, but these may vary across different datasets.
Warning 4: Do not run bash ipda_HyDRA_docker.sh blindly:
            - Check outputs of each step before proceeding.
            - Some steps (e.g., read processing) can run in parallel, e.g. short|long, but always verify step requirements before continuing.

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
"
}

#  ____       _   _   _                 
# / ___|  ___| |_| |_(_)_ __   __ _ ___ 
# \___ \ / _ \ __| __| | '_ \ / _` / __|
#  ___) |  __/ |_| |_| | | | | (_| \__ \
# |____/ \___|\__|\__|_|_| |_|\__, |___/
#                             |___/     

project_output_dir="path/from/working/dir/to/Project_Name/"
short_raw_dir="path/from/working/dir/to/data_short/" # Raw short_*1.f* and short_*2.f* files in fastq, fastq.gz, fq and fq.gz accepted
long_raw_dir="path/from/working/dir/to/data_long/" # Raw long_read.f* file(s) in fastq, fastq.gz, fq and fq.gz accepted

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
docker run --rm \
  -v "${short_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    fastqc -t 1 --outdir /output/hydra01S1_qc-raw_FastQC_${thislogdate} /data/short_read_1.fastq /data/short_read_2.fastq --memory 10000\
    "

## STEP 01S2 [resources:  -m 1 -c 1 -w "01:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    multiqc /output/hydra01S1_qc-raw_FastQC_DATE -o /output/hydra01S2_qc-raw_MultiQC_${thislogdate}\
    "

###############################################################
#                                                             #
#              STEP 01 [LONG-READS]: QC RAW FILES             #
#                                                             #
###############################################################

## STEP 01L1 [resources: -m 5 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${long_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    long_raw_stem=long-read-stem && \
    NanoPlot -t 1 -o /output/hydra01L1_qc-raw_NanoPlot_${thislogdate}/${long_raw_stem} --fastq /data/long-read-stem.fastq\
    "

## STEP 01L2 [resources: -m 5 -c 1 -w "04:00:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${long_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    names="stem1 stem2" && \
    NanoComp -t 1 --fastq $(ls /data/ | tr '\n' ' ') --names ${names} -o /output/hydra01L2_qc-raw_NanoComp_${thislogdate}
    "

## STEP 01L3 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${long_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    fastqc -t 1 --outdir /output/hydra01L3_qc-raw_FastQC_${thislogdate} /data/long-read-stem.fastq --memory 10000\
    "

## STEP 01L4 [resources: -m 1 -c 1 -w "00:30:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    multiqc /output/hydra01L3_qc-raw_FastQC_DATE -o /output/hydra01L4_qc-raw_MultiQC_${thislogdate}\
    "

## STEP 01L5 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${long_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    ls /data/ | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; seqtk seq -a ${file} > /output/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta; done  && \
    ls /output/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 /tools/fasta_splitter.py -f ${file} -n 20000 -p /output/hydra01L5_size-split_FastaSplitter_${thislogdate}/${main}
    "

###############################################################
#                                                             #
#        STEP 02 [SHORT-READS]: SHORT READS CORRECTION        #
#                                                             #
###############################################################

## STEP 02S1 [resources: -m 50 -c 15 -w "35:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
# WARNING 3: Set flag_merge to yes if you have samples sequenced from multiple lanes
docker run --rm \
  -v "${short_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    flag_merge="yes"\
    if [[ $flag_merge = "yes" ]]; then zcat -c /data/short_read_1.fastq | gzip -c >> /output/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R1.fastq.gz ; fi\
    if [[ $flag_merge = "yes" ]]; then zcat -c /data/short_read_2.fastq | gzip -c >> /output/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R2.fastq.gz ; fi\
    if [[ $flag_merge = "yes" ]]; then f1="/output/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R1.fastq.gz"; f2="/output/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R2.fastq.gz"; else f1=${short_raw_1}; f2=${short_raw_2}; fi; perl /tools/rcorrector/run_rcorrector.pl -t 15 -od /output/hydra02S1_short-cor_Rcorrector_${thislogdate} -1 ${f1} -2 ${f2}\
    bash /tools/bbmap/reformat.sh in=/output/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R1.cor.fq.gz in2=/output/hydra02S1_short-cor_Rcorrector_${thislogdate}/short_read_R2.cor.fq.gz\
    "

## STEP 02S2 [resources: -m 1 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${short_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    name="short_read_stem"\
    python /tools/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 /output/hydra02S1_short-cor_Rcorrector_DATE/short_read_R1.cor.fq.gz -2 /output/hydra02S1_short-cor_Rcorrector_DATE/short_read_R2.cor.fq.gz -s ${name}\
    mv /output/rmunfixable_${name}.log /output/hydra02S2_filter-uncor_FUPER_${thislogdate}/\
    mv /output/unfixrm_${name}*1.cor.fq /output/hydra02S2_filter-uncor_FUPER_${thislogdate}\
    mv /output/unfixrm_${name}*2.cor.fq /output/hydra02S2_filter-uncor_FUPER_${thislogdate}\
    "

###############################################################
#                                                             #
#          STEP 02 [LONG-READS]: LONG READS TRIMMING          #
#                                                             #
###############################################################

## STEP 02L1 [resources: -m 160 -c 40 -w "24:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${long_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    porechop -i /data/long-read-stem.fastq --threads 40 -o /output/hydra02L1_trim-adapters_Porechop_${thislogdate}/adapter-trimmed_long-read-stem.fastq.gz\
    "

## STEP 02L2 [resources: -m 2 -c 2 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    gunzip -c /output/hydra02L1_trim-adapters_Porechop_DATE/adapter-trimmed_long_raw_stem.fastq.gz | chopper -q 10 --headcrop 12 --threads 2 | gzip > /output/hydra02L2_trim-qc-head_Chopper_${thislogdate}/quality-trimmed_adapter-trimmed_long_raw_stem.fastq.gz\
    "

## STEP 02L3 [resources: -m 1 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    cutadapt --poly-a /output/hydra02L2_trim-qc-head_Chopper_DATE/quality-trimmed_adapter-trimmed_long_raw_stem.fastq.gz -o /output/hydra02L3_trim-polyA_Cutadapt_${thislogdate}/polyA-trimmed_quality-trimmed_adapter-trimmed_long_raw_stem.fastq.gz\
    "

###############################################################
#                                                             #
#          STEP 03 [SHORT-READS]: QC CORRECTED FILES          #
#                                                             #
###############################################################

## STEP 03S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    fastqc -t 1 --outdir /output/hydra03S1_qc-cor_FastQC_${thislogdate} /output/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_1.cor.fq /output/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_1.cor.fq --memory 10000\
    "

## STEP 03S2 [resources:  -m 1 -c 1 -w "01:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    multiqc /output/hydra03S1_qc-cor_FastQC_DATE -o /output/hydra03S2_qc-cor_MultiQC_${thislogdate}\
    "

###############################################################
#                                                             #
#           STEP 03 [LONG-READS]: QC TRIMMED FILES            #
#                                                             #
###############################################################

## STEP 03L1 [resources: -m 5 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    long_trim_stem=long_trim_stem && \
    NanoPlot -t 1 -o /output/hydra03L1_qc-trim_NanoPlot_${thislogdate}/${long_trim_stem} --fastq /output/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fastq.gz\
    "

## STEP 03L2 [resources: -m 5 -c 1 -w "10:00:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    names="stem1 stem2" && \
    NanoComp -t 1 --fastq $(ls /output/hydra02L3_trim-polyA_Cutadapt_DATE/ | tr '\n' ' ') --names ${names} -o /output/hydra03L2_qc-trim_NanoComp_${thislogdate}
    "

## STEP 03L3 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    fastqc -t 1 --outdir /output/hydra03L3_qc-trim_FastQC_${thislogdate} /output/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fastq.gz --memory 10000\
    "

## STEP 03L4 [resources: -m 1 -c 1 -w "00:30:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    multiqc /output/hydra03L3_qc-trim_FastQC_DATE -o /output/hydra03L4_qc-trim_MultiQC_${thislogdate}\
    "

## STEP 03L5 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING: Accepted formats: fastq, fastq.gz, fq and fq.gz
docker run --rm \
  -v "${long_raw_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    ls /data/hydra02L3_trim-polyA_Cutadapt_DATE/ | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; seqtk seq -a ${file} > /output/hydra03L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta; done  && \
    ls /output/hydra03L5_size-split_FastaSplitter_${thislogdate}/${main}.fasta | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 /tools/fasta_splitter.py -f ${file} -n 20000 -p /output/hydra03L5_size-split_FastaSplitter_${thislogdate}/${main}
    "

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
adapters_fasta_dir="/path/from/working/dir/to/adapters_trimmomatic_dir"
docker run --rm \
  -v "${adapters_fasta_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    trimmomatic PE -phred33 -threads 6 /output/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_stem_R1.cor.fq /output/hydra02S2_filter-uncor_FUPER_DATE/unfixrm_stem_R2.cor.fq /output/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/adapter-trimmed_unfixrm_stem_R1.fq.gz /output/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/unpaired_adapter-trimmed_unfixrm_stem_R1.fq.gz /output/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/adapter-trimmed_unfixrm_stem_R2.fq.gz /output/hydra04S1_trim-adapters_Trimmomatic_${thislogdate}/unpaired_adapter-trimmed_unfixrm_stem_R2.fq.gz ILLUMINACLIP:/data/adapters_trimmomatic.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:75 HEADCROP:12\
    "

## STEP 04S2 [resources: -m 5 -c 15 -w "00:30:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    /tools/bbmap/bbduk.sh t=15 -Xmx5g in=/output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R1.fq.gz in2=/output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R2.fq.gz qtrim=rl trimq=30 minavgquality=20 minbasequality=10 minlen=75 outm=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-failed_stem_R1.cor.fq.gz outm2=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-failed_stem_R2.cor.fq.gz bhist=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/stem_bhist.txt qhist=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_qhist.txt gchist=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_gchist.txt aqhist=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_aqhist.txt lhist=/output/hydra04S2_qc-ends-avg_BBDuk_${thislogdate}/quality-check_stem_lhist.txt gcbins=auto\
    "

###############################################################
#                                                             #
#         STEP 04 [LONG-READS]: LONG READS CORRECTION         #
#                                                             #
###############################################################

## STEP 04L1 [resources: -m 60 -c 5 -w "10:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    gunzip -c /output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_* | awk 'NR % 4 == 2' | tr NT TN | ropebwt2 -LR | tr NT TN >> /output/hydra04L1_short-reads-bwt_RopeBWT2_${thislogdate}/raw_all-short-reads.bwt\
    "

## STEP 04L2 [resources: -m 20 -c 2 -w "25:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    fmlrc2-convert -i /output/hydra04L1_short-reads-bwt_RopeBWT2_DATE/raw_all-short-reads.bwt /output/hydra04L2_short-reads-bwt_FMLRC2-convert_${thislogdate}/comp_msbwt_all-short-reads.npy\
    "

## STEP 04L3 [resources: -m 300 -c 8 -w "12:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    fmlrc2 --threads 8 --cache_size 8 /output/hydra04L2_short-reads-bwt_FMLRC2-convert_DATE/comp_msbwt_all-short-reads.npy /output/hydra02L3_trim-polyA_Cutadapt_DATE/polyA-trimmed_quality-trimmed_adapter-trimmed_long-read-stem.fastq.gz /output/hydra04L3_long-cor_FMLRC2_${thislogdate}/corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_long-read-stem.fa\
    gzip /output/hydra04L3_long-cor_FMLRC2_${thislogdate}/corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_long-read-stem.fa\
    "

###############################################################
#                                                             #
#           STEP 05 [SHORT-READS]: QC TRIMMED FILES           #
#                                                             #
###############################################################

## STEP 05S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    fastqc -t 1 --outdir /output/hydra05S1_qc-trim_FastQC_${thislogdate} /output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_1.fq.gz /output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_2.fq.gz --memory 10000\
    "

## STEP 05S2 [resources:  -m 1 -c 1 -w "01:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    multiqc /output/hydra05S1_qc-trim_FastQC_DATE -o /output/hydra05S2_qc-trim_MultiQC_${thislogdate}\
    "

###############################################################
#                                                             #
#           STEP 05 [LONG-READS]: QC CORRECTED FILES          #
#                                                             #
###############################################################

## STEP 05L1 [resources: -m 1 -c 1 -w "00:30:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    ls /data/hydra04L3_long-cor_FMLRC2_DATE/ | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; gunzip -c ${file}.fa.gz > /output/hydra05L1_size-split_FastaSplitter_${thislogdate}/${main}.fasta; done  && \
    ls /output/hydra05L1_size-split_FastaSplitter_${thislogdate}/${main}.fasta | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 /tools/fasta_splitter.py -f ${file} -n 20000 -p /output/hydra05L1_size-split_FastaSplitter_${thislogdate}/${main}
    "

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
ribosomal_rna_ref_dir="/path/from/working/dir/to/ribosomalRNA_dir"
docker run --rm \
  -v "${ribosomal_rna_ref_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    /tools/bbmap/bbduk.sh t=6 -Xmx1g in=/output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R1.fq.gz in2=/output/hydra04S1_trim-adapters_Trimmomatic_DATE/adapter-trimmed_unfixrm_stem_R2.fq.gz outm=/output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribosomal-RNA_stem_R1.cor.fq.gz outm2=/output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribosomal-RNA_stem_R2.cor.fq.gz out=/output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R1.cor.fq out2=/output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R2.cor.fq ref=/data/ribosomalRNA.fa stats=/output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/rRNA-alignment_stats.txt\
    cat /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R1.cor.fq | awk '{ if (NR%4==1) { print $1substr($2,2)"/1" } else { print } }' | sed 's/_1://g' >> /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R1.cor.fq; gzip /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R1.cor.fq; gzip /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R1.cor.fq\
    cat /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R2.cor.fq | awk '{ if (NR%4==1) { print $1substr($2,2)"/2" } else { print } }' | sed 's/_2://g' >> /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R2.cor.fq; gzip /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/header-fixed_ribo-depleted_stem_R2.cor.fq; gzip /output/hydra06S1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_stem_R2.cor.fq\
    "

###############################################################
#                                                             #
#           STEP 06 [LONG-READS]: FILTER rRNA READS           #
#                                                             #
###############################################################

## STEP 06L1 [resources: -m 5 -c 5 -w "00:30:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ribosomal RNAs file
# WARNING 3: BBDuk uses more memory and CPUs than what is required. Therefore, we only ask for 60% of the memory and CPU specified when running this step.
ribosomal_rna_ref_dir="/path/from/working/dir/to/ribosomalRNA_dir"
docker run --rm \
  -v "${ribosomal_rna_ref_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    /tools/bbmap/bbduk.sh t=3 -Xmx3g in=/output/hydra04L3_long-cor_FMLRC2_DATE/corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fa.gz outm=/output/hydra06L1_ribodepletion_BBDuk_${thislogdate}/ribosomal-RNA_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fasta out=/output/hydra06L1_ribodepletion_BBDuk_${thislogdate}/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fasta ref=/data/ribosomalRNA.fa stats=/output/hydra06L1_ribodepletion_BBDuk_${thislogdate}/rRNA-alignment_stats.txt\
    "

###############################################################
#                                                             #
#          STEP 07 [SHORT-READS]: MAP READS 2 GENOME          #
#                                                             #
###############################################################

## STEP 07S1 [resources: -m 30 -c 15 -w "02:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to bowtie2 ref genome index
genome_bowtie2_index_dir="path/from/working/dir/to/ref_genome_bowtie2_dir"
docker run --rm \
  -v "${genome_bowtie2_index_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    bowtie2 --fast --no-unal -p 15 -x /data/ref_genome_bowtie2 -1 /output/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R1.cor.fq.gz -2 /output/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R2.cor.fq.gz | samtools view -@15 -bS - >> /output/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.bam\
    samtools sort -@15 /output/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.bam -o /output/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam\
    samtools index -@15 /output/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam /output/hydra07S1_gen-map_Bowtie2_${thislogdate}/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam.bai\
    "

###############################################################
#                                                             #
#          STEP 07 [LONG-READS]: MAP READS 2 GENOME           #
#                                                             #
###############################################################

## STEP 07L1 [resources: -m 150 -c 15 -w "30:00:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to bowtie2 ref genome index
genome_bowtie2_index_dir="path/from/working/dir/to/ref_genome_bowtie2_dir"
docker run --rm \
  -v "${genome_bowtie2_index_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    bowtie2 --fast --no-unal -p15 -x /data/ref_genome_bowtie2 -f /output/hydra06L1_ribodepletion_BBDuk_DATE/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.fasta | samtools view -@15 -bS - >> /output/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.bam\
    samtools sort -@15 /output/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.bam -o /output/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.coordsorted.bam\
    samtools index -@15 /output/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.bam /output/hydra07L1_gen-map_Bowtie2_${thislogdate}/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.coordsorted.bam.bai\
    "

###############################################################
#                                                             #
#          STEP 08 [SHORT-READS]: STRAND ASSESSMENT           #
#                                                             #
###############################################################

## STEP 08S1 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ref gene BED file. e.g. download <https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz>, unzip files and use sed to change chromosome nomenclature (to match your reference genome): sed 's/^chr//g'. WARNING! save files as hg38RefGene.bed
reference_gene_dir="/path/from/working/dir/to/ref_gene_dir"
docker run --rm \
  -v "${reference_gene_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    infer_experiment.py -r /data/ref_gene.bed -i /output/hydra07S1_gen-map_Bowtie2_DATE/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.coordsorted.bam >> /output/hydra08S1_strandness-assess_RSeQC_${thislogdate}/strand-assessment_genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem.txt\
    "

## STEP 08S2 [resources: -m 1 -c 8 -w "02:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    hex_name_colors='black;grey;lightgrey;palevioletred;firebrick;crimson;orangered;saddlebrown;bisque;navajowhite;orchid;magenta;darkkhaki;olive;darkolivegreen;darkseagreen;darkgreen;mediumseagreen;aquamarine;lightcyan;darkcyan;cadetblue;lightskyblue;lightslategrey;royalblue;darkblue;darkslateblue;indigo;plum;fuchsia;hotpink;lightpink'\
    colors=`echo ${hex_name_colors} | tr ';' '\n' | tr '\n' ' ' | head -n${sample_num}`\
    labels="stem1 stem2"\
    echo '${hex_name_colors}' | tr ';' '\n' >> /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-pallete_${thislogdate}.txt\
    echo '${labels}' | tr ' ' '\n' | uniq >> /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels_${thislogdate}.txt\
    head -n${sample_num} /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-pallete_${thislogdate}.txt >> /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-scheme_${thislogdate}.txt\
    paste /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels_${thislogdate}.txt /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-color-scheme_${thislogdate}.txt >> /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels-color-scheme_${thislogdate}.txt\
    cut -f2 /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-labels-color-scheme_${thislogdate}.txt >> /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-sample-colors_${thislogdate}.txt\
    multiBamSummary bins --bamfiles /output/hydra07S1_gen-map_Bowtie2_DATE/genome-alignment_header-fixed_ribo-depleted_adapter-trimmed_unfixrm*.coordsorted.bam --labels ${labels} --outFileName /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --maxFragmentLength 280 --numberOfProcessors 8\
    plotCorrelation --corData /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --corMethod spearman --whatToPlot heatmap --labels ${labels} --plotFile /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/spearman-correlation_analysis-stem.pdf --plotTitle Spearman-Correlation_analysis-stem_${thislogdate} --outFileCorMatrix /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/spearman-correlation_analysis-stem.matrix --plotNumbers\
    plotCorrelation --corData /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --corMethod pearson --whatToPlot heatmap --labels ${labels} --plotFile /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pearson-correlation_analysis-stem.pdf --plotTitle Pearson-Correlation_analysis-stem_${thislogdate} --outFileCorMatrix /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pearson-correlation_analysis-stem.matrix --plotNumbers\
    plotPCA --corData /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/bam-comparison-summary_analysis-stem.npz --plotHeight 20 --labels ${labels} --colors ${colors} --plotFile /output/hydra08S2_rep-cor_DeepTools_${thislogdate}/pca-plot_analysis-stem.pdf --plotTitle PCA-plot_analysis-stem_${thislogdate} --log2\
    "

###############################################################
#                                                             #
#           STEP 08 [LONG-READS]: STRAND ASSESSMENT           #
#                                                             #
###############################################################

## STEP 08L1 [resources: -m 1 -c 1 -w "00:30:00"]
# WARNING 1: You will need to loop through sample files.
# WARNING 2: User must provide path to ref gene BED file. e.g. download <https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz>, unzip files and use sed to change chromosome nomenclature (to match your reference genome): sed 's/^chr//g'. WARNING! save files as hg38RefGene.bed
reference_gene_dir="/path/from/working/dir/to/ref_gene_dir"
docker run --rm \
  -v "${reference_gene_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    infer_experiment.py -r /data/ref_gene.bed -i /output/hydra07L1_gen-map_Bowtie2_DATE/genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.coordsorted.bam >> /output/hydra08L1_strandness-assess_${thislogdate}/strand-assessment_genome-alignment_ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_stem.txt\
    "

###############################################################
#                                                             #
#        STEP 09 [SHORT-READS]: QC RIBODEPLETED FILES         #
#                                                             #
###############################################################

## STEP 09S1 [resources: -m 10 -c 1 -w "02:00:00"]
# WARNING: You will need to loop through sample files.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    fastqc -t 1 --outdir /output/hydra09S1_qc-ribodepleted_FastQC_${thislogdate} /output/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_1.cor.fq.gz /output/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_2.cor.fq.gz --memory 10000\
    "

## STEP 09S2 [resources:  -m 1 -c 1 -w "01:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \ 
    multiqc /output/hydra09S1_qc-ribodepleted_FastQC_DATE -o /output/hydra09S2_qc-ribodepleted_MultiQC_${thislogdate}\
    "

###############################################################
#                                                             #
#         STEP 09 [LONG-READS]: QC RIBODEPLETED FILES         #
#                                                             #
###############################################################

## STEP 09L1 [resources: -m 1 -c 1 -w "00:30:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    ls /output/hydra06L1_ribodepletion_BBDuk_DATE/${main}.fasta | while read file; do main=`echo "$(basename "${file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; python3 /tools/fasta_splitter.py -f ${file} -n 20000 -p /output/hydra09L1_size-split_FastaSplitter_${thislogdate}/${main}\
    "

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
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    zcat -c /output/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm*_R1.cor.fq.gz | gzip -c >> /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R1.cor.fq.gz\
    zcat -c /output/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm*_R2.cor.fq.gz | gzip -c /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R2.cor.fq.gz\
    gzip -c /output/hydra06L1_ribodepletion_BBDuk_DATE/ribo-depleted_corrected_polyA-trimmed_quality-trimmed_adapter-trimmed_*.fasta >> /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-long-reads_merged.fa.gz\
    rnaspades.py --checkpoint "all" -t 40 -m 250 -1 /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R1.cor.fq.gz -2 /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R2.cor.fq.gz --nanopore /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-long-reads_merged.fa.gz -o /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}\
    seqtk seq /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/transcripts.fasta > /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/transcripts.singleline.fasta\
    rm -f /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R1.cor.fq.gz /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-short-reads_merged-paired-end_R2.cor.fq.gz /output/hydra11H1_assembly_rnaSPAdes_${thislogdate}/all-long-reads_merged.fa.gz\
    "

## STEP 11H2 [resources: -m 250 -c 12 -w "85:00:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    ls ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R1.cor.fq.gz >> trinity-left-list_${thislogdate}.txt\
    ls ${project_output_dir}/hydra06S1_ribodepletion_BBDuk_DATE/header-fixed_ribo-depleted_adapter-trimmed_unfixrm_stem_R2.cor.fq.gz >> trinity-right-list_${thislogdate}.txt\
    left_list=`cat trinity-left-list_${thislogdate}.txt | tr '\n' ',' | tr -d ' '`; left_list="${left_list%\,}"; right_list=`cat trinity-right-list_${thislogdate}.txt | tr '\n' ',' | tr -d ' '`; right_list="${right_list%\,}"\
    Trinity --seqType fq --max_memory 250G --left $left_list --right $right_list --CPU 12 --SS_lib_type RF --output ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate} --group_pairs_distance 280 --full_cleanup --bflyHeapSpaceMax 16G --normalize_max_read_cov 50 --bypass_java_version_check\
    mv trinity-left-list_${thislogdate}.txt trinity-right-list_${thislogdate}.txt ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}\
    mv ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}.Trinity.fasta  ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}/assembly.Trinity.fasta\
    mv ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}.Trinity.fasta.gene_trans_map ${project_output_dir}/hydra11H2_assembly_Trinity_${thislogdate}/assembly.Trinity.fasta.gene_trans_map\
    "

###############################################################
#                                                             #
#          STEP 12 [HYBRID]: QUALITY CHECK ASSEMBLY           #
#                                                             #
###############################################################

## STEP 12H1 [resources: -m 2 -c 5 -w "50:00:00"]
# WARNING: User must provide path to BUSCO lineage dataset, which can be obtained from <https://busco-archive.ezlab.org/v3/datasets/eukaryota_odb9.tar.gz>. Filename eukaryota_odb9 given as an example
busco_odb9_dir="/path/from/working/dir/to/BuscoLineages/"
docker run --rm \
  -v "${busco_odb9_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    BUSCO.py -i /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta -l /data/eukaryota_odb9 -o eukaryota-odb9_stem-rnaSPAdes -m transcriptome -c 5 -t /output/hydra12H1_qc_BUSCO_${thislogdate} -sp human\
    BUSCO.py -i /output/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta -l /data/eukaryota_odb9 -o eukaryota-odb9_stem-Trinity -m transcriptome -c 5 -t /output/hydra12H1_qc_BUSCO_${thislogdate} -sp human\
    "

## STEP 12H2 [resources: -m 1 -c 1 -w "00:30:00"]
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    ${trinstats} /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta >> /output/hydra12H2_qc_TrinityStats_${thislogdate}/trinityStats_stem-rnaSPAdes\
    ${trinstats} /output/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta >> /output/hydra12H2_qc_TrinityStats_${thislogdate}/trinityStats_stem-Trinity\
    "

## STEP 12H3 [resources: -m 80 -c 15 -w "120:00:00"]
# WARNING: User must provide path to reference transcriptome annotation file
reference_fasta_dir="/path/from/working/dir/to/ref-transcriptome-dir/"
docker run --rm \
  -v "${reference_fasta_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    transrate --assembly /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta --reference /data/ref-transcriptome.fasta --threads 15 --output /output/hydra12H3_qc_TransRate_${thislogdate}/\
    mv /output/hydra12H3_qc_TransRate_${thislogdate}/assemblies.csv /output/hydra12H3_qc_TransRate_${thislogdate}/stem_assemblies.csv\
    "

## STEP 12H4 [resources: -m 85 -c 32 -w "05:00:00"]
# WARNING: User must provide path to reference transcriptome annotation file
reference_fasta_dir="/path/from/working/dir/to/ref-transcriptome-dir/"
docker run --rm \
  -v "${reference_fasta_dir}":/data
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    tablesize=`awk 'END{print NR}' hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq | cut -d' ' -f1`\
    fastq_pair -t ${tablesize} hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/right.norm.fq\
    transrate --assembly /output/hydra11H1_assembly_rnaSPAdes_DATE/assembly.Trinity.fasta --left hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq --right hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/right.norm.fq.paired.fq --reference ${reference_fasta} --threads 32 --output /output/hydra12H4_qc-and-polish_TransRate_${thislogdate}\
    mv /output/hydra12H4_qc-and-polish_TransRate_${thislogdate}/assemblies.csv /output/hydra12H4_qc-and-polish_TransRate_${thislogdate}/stem_assemblies.csv\
    "

###############################################################
#                                                             #
#     STEP 13 [HYBRID]: TRANSCRIPTOME READ REPRESENTATION     #
#                                                             #
###############################################################

## STEP 13H1 [resources: -m 250 -c 16 -w "120:00:00"]
# WARNING: You will need to loop through sample files for Bowtie2 and GMAP alignments.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    bowtie2-build --threads 16 /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_bowtie2idx\
    bowtie2 -p 16 -q --no-unal -k 20 -x /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_bowtie2idx -1 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R1.cor.fq.gz -2 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R2.cor.fq.gz | samtools view -@16 -Sb -o /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.bam\
    gmap -n 0 -t 16 -B 5 --dir=/output/hydra11H1_assembly_rnaSPAdes_DATE/ --gseg /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta /output/hydra06L1_ribodepletion_BBDuk_DATE/stem.fasta | samtools view -@16 -bS - >> /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.bam\
    samtools merge -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.bam /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.bam /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.bam\
    samtools sort -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.bam -o /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.coordsorted.bam\
    samtools sort -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.bam -o /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam\
    samtools sort -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.bam -o /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.coordsorted.bam\
    samtools index -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.coordsorted.bam\
    samtools index -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam\
    samtools index -@16 /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.coordsorted.bam\
    samtools faidx /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta -o /output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_samtoolsidx\
    ${picardtool} CollectInsertSizeMetrics I=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_merged_shortreads-bowtie2-and-longreads-gmap.coordsorted.bam O=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-transcriptome-alignment-merged.txt H=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-merged-histogram.pdf INCLUDE_DUPLICATES=FALSE\
    ${picardtool} CollectInsertSizeMetrics I=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam O=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-transcriptome-alignment-shortreads.txt H=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-shortreads-histogram.pdf INCLUDE_DUPLICATES=FALSE\
    ${picardtool} CollectInsertSizeMetrics I=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_transcriptome-alignment-longreads-gmap.coordsorted.bam O=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-transcriptome-alignment-longreads.txt H=/output/hydra13H1_read-representation_GMAP-Bowtie2_${thislogdate}/transcripts_insert-size-longreads-histogram.pdf INCLUDE_DUPLICATES=FALSE\
    "

## STEP 13H1 [resources: -m 250 -c 16 -w "120:00:00"]
# WARNING: You will need to loop through sample files for Bowtie2 alignments.
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    bowtie2-build --threads 15 /output//output/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_bowtie2idx\
    bowtie2 -p 15 -q --no-unal -k 20 -x /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_bowtie2idx -1 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R1.cor.fq.gz -2 HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R2.cor.fq.gz | samtools view -@15 -Sb -o /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.bam\
    samtools sort -@15 /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.bam -o /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam\
    samtools index -@15 /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-shortreads-bowtie2.coordsorted.bam\
    samtools faidx /output//output/hydra11H2_assembly_Trinity_DATE/assembly.Trinity.fasta -o /output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_samtoolsidx\
    ${picardtool} CollectInsertSizeMetrics I=/output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_transcriptome-alignment-longreads-gmap.coordsorted.bam O=/output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_insert-size-transcriptome-alignment-longreads.txt H=/output/hydra13H2_read-representation_Bowtie2_${thislogdate}/assembly_insert-size-longreads-histogram.pdf INCLUDE_DUPLICATES=FALSE\
    "

###############################################################
#                                                             #
#            STEP 14 [HYBRID]: SPLICING ASSESSMENT            #
#                                                             #
###############################################################

## STEP 14H1 [resources: -m 15 -c 16 -w "10:00:00"]
# WARNING 1: User must provide name of reference_GMAP_genome_db_dir and its path separately.
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
thislogdate=$(date +'%d%m%Y%H%M%S%Z') 
referencedir=/path/from/working/dir/to/
docker run --rm \
  -v "${referencedir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    gmap -n 0 -t 16 -B 5 --dir=/data --db=reference_GMAP_genome_db_dir /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts_splicing-assessment.genome.sam\
    grep -v "^@" /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts_splicing-assessment.genome.sam | awk -v A="/output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.duplicated-alignments.sam" -v B="/output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam" '{count[$1]++; lines[$1,count[$1]] = $0;} END {for (key in count) {file = (count[key] == 2 ? A : B); for (i = 1; i <= count[key]; i++) print lines[key, i] > file;}}'\
    sed -n 'p;n' /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.duplicated-alignments.sam >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.sam\
    cat /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam | awk '$1 ~ /@/ || ($6 !~ /N/ && $6 !~ /\*/)' >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.sam\
    cat /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam | awk '$1 ~ /@/ || ($6 ~ /N/)' >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.sam\
    cat /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.single-alignments.sam | awk '$1 ~ /@/ || ($6 ~ /\*/)' >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unaligned.sam\
    grep -v "^@" /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.sam | cut -f1 >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.ids\
    grep -v "^@" /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.sam | cut -f1 >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.ids\
    seqtk subseq /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.ids >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.spliced.fasta\
    seqtk subseq /output/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.fasta /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.ids >> /output/hydra14H1_splicing-assessment_GMAP_${thislogdate}/transcripts/transcripts.unspliced.fasta\
    "

###############################################################
#                                                             #
#                STEP 15 [HYBRID]: READ SUPORT                #
#                                                             #
###############################################################

## STEP 15H1 [resources: -m 20 -c 16 -w "20:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: RPKM cutoffs set to 1 (high) and 0.5 (low) below. Change accordingly
# WARNING 3: CDhit cutoffs set to coverage of 90 and identity of 98. Change accordingly
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    high_rpkm_cutoff=1\
    low_rpkm_cutoff=0.5\
    cd-hit-est -o /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -c 0.98 -i /output/hydra14H1_splicing-assessment_GMAP_DATE/transcripts/transcripts.spliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.90 -p 1 -d 0 -b 3 -T 16 -M 20000\
    minimap2 -ax splice -t 16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra06L1_ribodepletion_BBDuk_DATE/stem1.fasta,/output/hydra06L1_ribodepletion_BBDuk_DATE/stem2.fasta | samtools view -@16 -bS - >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.bam\
    bowtie2-build --threads 16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx\
    bowtie2 -q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --threads 16 -x /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx -1 /output/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq  -2 /output/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq  | samtools view -@16 -bS - >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam\
    samtools merge /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap.spliced.bam -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.bam\
    samtools sort -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam -o /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.coordsorted.bam\
    samtools sort -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.bam -o /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.coordsorted.bam\
    samtools sort -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap.spliced.bam -o /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap.spliced.coordsorted.bam\
    rm -f /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.bam /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap.spliced.bam\
    samtools index -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.coordsorted.bam\
    samtools index -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.coordsorted.bam\
    samtools index -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap.spliced.coordsorted.bam\
    samtools idxstats -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.coordsorted.bam >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts.spliced.tsv\
    samtools idxstats -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap.spliced.coordsorted.bam >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap_alignment-counts.spliced.tsv\
    samtools idxstats -@16 /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap.spliced.coordsorted.bam >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap_alignment-counts.spliced.tsv\
    total=`sed '$d' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv\
    total=`sed '$d' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv\
    total=`sed '$d' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv\
    awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta\
    awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta\
    awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta\
    awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta\
    awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_longreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta\
    awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_shortreads_bt2_and_longreads_minimap_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> /output/hydra15H1_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_merged_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta\
    "

## STEP 15H2 [resources: -m 20 -c 16 -w "20:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: RPKM cutoffs set to 1 (high) and 0.5 (low) below. Change accordingly
# WARNING 3: CDhit cutoffs set to coverage of 90 and identity of 98. Change accordingly
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    high_rpkm_cutoff=1\
    low_rpkm_cutoff=0.5\
    cd-hit-est -o /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -c 0.98 -i /output/hydra14H1_splicing-assessment_GMAP_DATE/transcripts/transcripts.spliced.fasta -g 1 -r 0 -aL 0.9 -aS 0.90 -p 1 -d 0 -b 3 -T 16 -M 20000\
    bowtie2-build --threads 16 /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx\
    bowtie2 -q --no-mixed --no-discordant --gbar 1000 --np 0 --end-to-end -k 200 --threads 16 -x /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_spliced_bt2idx -1 /output/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq -2 /output/hydra11H2_assembly_Trinity_DATE/insilico_read_normalization/left.norm.fq.paired.fq | samtools view -@16 -bS - >> /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam\
    samtools sort -@16 /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam -o /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.coordsorted.bam\
    rm -f /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.bam\
    samtools index -@16 /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.coordsorted.bam\
    samtools idxstats -@16 /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2.spliced.coordsorted.bam >> /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts.spliced.tsv\
    total=`sed '$d' /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts.spliced.tsv | awk 'BEGIN {OFS="\t"} {total_reads += $3} END {print total_reads}'` ; sed '$d' /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts.spliced.tsv | awk -v total_reads="${total}" 'BEGIN {OFS="\t"} {scaling_factor = total_reads / 1000000} {rpm = ($3 == 0) ? 0 : $3 / scaling_factor} {rpkm = ($2 == 0) ? 0 : (rpm == 0) ? 0 : (rpm / ($2 / 1000))} {print $1, $2, $3, $4, rpm, rpkm}' > /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv\
    awk -v rpkm_cutoff="${high_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.ids >> /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_highcutoff${high_rpkm_cutoff}.spliced.fasta\
    awk -v rpkm_cutoff="${low_rpkm_cutoff}" 'NR>0{if ($6 > rpkm_cutoff) print}' /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_bt2_alignment-counts_rpm_rpkm.spliced.tsv | cut -f1 >> /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids ; seqtk subseq /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.ids >> /output/hydra15H2_read-support_Samtools_${thislogdate}/transcripts_spliced/transcripts_CDhit_identity98_cov90_shortreads_RPKM_lowcutoff${low_rpkm_cutoff}.spliced.fasta\
    "

###############################################################
#                                                             #
#         STEP 16 [HYBRID]: TRANSCRIPTOME ANNOTATION          #
#                                                             #
###############################################################

## STEP 16H1 [resources: -m 20 -c 16 -w "10:00:00"]
# WARNING 1: User must provide: name of reference_GMAP_genome_db_dir and its path separately; reference transcriptome BED file; lncRNA database BED file
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
# WARNING 3: Please set
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
referencedir=/path/from/working/dir/to/
transcriptome_bed_dir="/path/from/working/dir/to/ref_transcriptome_dir"
lncrnadb_dir="/path/from/working/dir/to/lncRNAdb_dir"
docker run --rm \
  -v "${referencedir}":/data_gmap \
  -v "${transcriptome_bed_dir}":/data_transcriptome \
  -v "${lncrnadb_dir}":/data_lncRNA \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    gmap -n 1 -t 16 -B 5 --dir=/data_gmap --db=reference_GMAP_genome_db_dir /output/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta --min-intronlength=50 --max-intronlength-middle=200000 --min-identity=0.9 --format=samse > /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.sam\
    convert2bed-megarow --input=sam --do-not-sort --max-mem=20G < /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.sam > /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.bed12\
    cut -f1-6 /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.bed12 >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.bed\
    bedtools intersect -wo -s -f 0.75 -r -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.bed -b /data_transcriptome/ref-transcriptome.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomefeatures.bed\
    bedtools intersect -s -f 0.75 -r -C -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.bed -b /data_transcriptome/ref-transcriptome.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomecounts.bed\
    bedtools intersect -v -wa -s -f 0.75 -r -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.bed -b /data_transcriptome/ref-transcriptome.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.notingenome.bed\
    bedtools intersect -wo -s -f 0.75 -r -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.notingenome.bed -b /data_lncRNA/lncRNAdb.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.lncRNAfeatures.bed\
    bedtools intersect -s -f 0.75 -r -C -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.notingenome.bed -b /data_lncRNA/lncRNAdb.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.lncRNAcounts.bed\
    bedtools intersect -v -wa -s -f 0.75 -r -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.notingenome.bed -b /data_lncRNA/lncRNAdb.bed /data_transcriptome/ref-transcriptome.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.novel.bed\
    bedtools intersect  -wo -s -f 0.5 -e -a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.novel.bed -b /data_lncRNA/lncRNAdb.bed /data_transcriptome/ref-transcriptome.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.potentialfragment.bed\
    cut -f4 /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomefeatures.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_a ; cut -d'"' -f2,6 /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomefeatures.bed | tr '"' '\t' >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_b ; grep -Po 'gene_biotype "\K.*?(?=")' /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomefeatures.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_c ; grep -Po 'transcript_biotype "\K.*?(?=")' /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomefeatures.bed >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_d ; paste /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_b /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_c /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_d >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation; rm -rf /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_a /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_b /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_c /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.tmp_d ; cut -f4,10,11 /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.lncRNAfeatures.bed | sed "s/\t./\.\tncRNAdb\t\./g" >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation\
    cut -f1,4 /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation | sort | uniq >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary\
    grep "protein_coding" /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary.protein.transcriptID\
    grep "lincRNA\|ncRNAdb\|lncRNA" /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary.ncRNA.transcriptID\
    grep "antisense" /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary.antisense.transcriptID\
    grep "pseudogene" /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.annotation.summary.pseudogene.transcriptID\
    cut -f12,16 /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem.genomefeatures.bed | tr '\"' '\t' | sort | uniq >> /output/hydra16H1_transcriptome-annotation_BedTools_${thislogdate}/transcripts/transcripts.refgenome_stem_geneID-type_retrieved-annotated-genes\
    "

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
genome_proteins_dir="/path/from/working/dir/to/protein-coding_reference-annotation_dir"
docker run --rm \
  -v "${genome_proteins_dir}":/data \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    bedtools intersect -wo -s -f 0.6 -r -a /output/hydra16H1_transcriptome-annotation_BedTools_DATE/transcripts/transcripts.bed -b /data/protein-coding_reference-annotation.bed >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.proteincoding.bed\
    ezLncPred -i /output/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -o /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2 CPC2\
    ezLncPred -i /output/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -o /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT CPAT -p Human\
    ezLncPred -i /output/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta -o /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI CNCI -p ve\
    cut -f1,7,8 /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2 | grep -w "coding" | sort | uniq >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.proteincoding\
    cut -f1,7,8 /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2 | grep -w "noncoding" | sort | uniq >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.noncoding\
    cat /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT | tr "'" '\t' | cut -f2,8 | awk 'NR>0{if($2>=0.5){print}}' | sort | uniq >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.proteincoding\
    cat /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT | tr "'" '\t' | cut -f2,8 | awk 'NR>0{if($2<0.5){print}}' | sort | uniq >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.noncoding\
    cat /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI/CNCI.index | tail -n +2 | cut -f1,2,3 | grep -Pe "\tcoding" | awk 'NR>0{printf "%s\t%.2f\n", $1,$3}' | sort | uniq >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI.proteincoding\
    cat /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI/CNCI.index | tail -n +2 | cut -f1,2,3 | grep -Pe "\tnoncoding" | awk 'NR>0{printf "%s\t%.2f\n", $1,$3}' | sort | uniq >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CNCI.noncoding\
    number=`cut -f1 ${annotation} /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.proteincoding | tail -n +2 | sort | uniq -c | grep -c " 2 "` ; echo "/output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts\tCPC2 Annotated ptcoding classified as ptcoding\t${number}" >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPC2.proteincodingtest\
    number=`cut -f1 ${annotation} /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.proteincoding | tail -n +2 | sort | uniq -c | grep -c " 2 "` ; echo "/output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts\tCPAT Annotated ptcoding classified as ptcoding\t${number}" >> /output/hydra18H1_coding-potential_ezLncPred_${thislogdate}/transcripts/transcripts.CPAT.proteincodingtest\
    "

###############################################################
#                                                             #
#         STEP 19 [HYBRID]: IDENTIFY lncRNA POTENTIAL         #
#                                                             #
###############################################################

## STEP 19H1 [resources: -m 120 -c 16 -w "60:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
# WARNING 3: User must provide: 
gencodegenome_dir="/path/from/working/dir/to/gencode.genome_dir" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz>
gencodebed_dir="/path/from/working/dir/to/gencode.v36.annotation_dir"
transgtf_dir="/path/from/working/dir/to/ref-transcriptome_dir"
reference_fasta_dir="/path/from/working/dir/to/ref-transcriptome.fasta"
ptngtf_dir="/path/from/working/dir/to/gencode.v36.proteincoding_dir" # Protein-coding transcripts were retrieved using grep "^#\|protein_coding" from the above file
lncrnagtf_dir="/path/from/working/dir/to/gencode.v36.long_noncoding_RNAs.gtf" # e.g. <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.long_noncoding_RNAs.gtf.gz>
goodlncrnagtf_dir="/path/from/working/dir/to/gencode.v36.confirmed_long_noncoding_RNAs_dir" # grep -v "TEC" out of above (to be experimentally confirmed transcripts)
docker run --rm \
  -v "${gencodegenome_dir}":/data_gencodegenome_dir \
  -v "${gencodebed_dir}":/data_gencodebed_dir \
  -v "${transgtf_dir}":/data_transgtf_dir \
  -v "${ptngtf_dir}":/data_ptngtf_dir \
  -v "${goodlncrnagtf_dir}":/data_goodlncrnagtf_dir \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    minimap2 -x splice --secondary=no -a -t 16 -o /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sam /data_gencodegenome_dir/gencode.genome.fa /output/hydra15H1_read-support_Samtools_DATE/transcripts_spliced/transcripts_CDhit_identity98_cov90.spliced.fasta\
    samtools view -u /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sam | samtools sort -@ 16 -o /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sorted.bam\
    bamToBed -splitD -bed12 -i /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.sorted.bam > /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.bed12\
    bedtools intersect -wo -s -f 0.75 -r -a /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.bed12 -b /data_gencodebed_dir/gencode.v36.annotation.bed >> /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.gencode.genomefeatures.bed\
    bedToGenePred /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.bed12 stdout | genePredToGtf file stdin /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.gtf\
    FEELnc_filter.pl -i /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.gtf -a /data_transgtf_dir/ref-transcriptome.gtf --biotype transcript_biotype=protein_coding --size=200 --monoex=0 -p 16 -o /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.feelncfilter.log > /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.lncRNAcandidates.gtf\
    FEELnc_codpot.pl -i /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/transcripts.lncRNAcandidates.gtf -a /data_ptngtf_dir/gencode.v36.proteincoding.gtf -l /data_goodlncrnagtf_dir/gencode.v36.confirmed_long_noncoding_RNAs.gtf -g /data_gencodegenome_dir/gencode.genome.fa -o transcripts --outdir=/output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts\
    FEELnc_classifier.pl -i /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/transcripts.lncRNA.gtf -a /data_ptngtf_dir/gencode.v36.proteincoding.gtf -l /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/transcripts.classifier.log >> /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/transcripts.lncRNAclasses.txt\
    mkdir -p /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/VarImpPlots /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/ROCcurves\
    mv /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/*varImpPlot.png /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/VarImpPlots ; mv /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts/codpot_transcripts/*TGROC.png /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/ROCcurves\
    "

###############################################################
#                                                             #
#              STEP 20 [HYBRID]: DEFINE lncRNAs               #
#                                                             #
###############################################################

## STEP 20H1 [resources: -m 1 -c 1 -w "01:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    grep -P "^1\t" /output/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt | cut -f3- >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.noncoding\
    cut -f1 /output/hydra18H1_coding-potential_ezLncPred_DATE/stem/stem.CNCI.noncoding /output/hydra18H1_coding-potential_ezLncPred_DATE/stem/stem.CPAT.noncoding /output/hydra18H1_coding-potential_ezLncPred_DATE/stem/stem.CPC2.noncoding | sort | uniq -i -c | grep -v " 1 " | cut -c9- >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred.min2-models.noncoding\
    cut -f1 /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred.min2-models.noncoding /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.noncoding | sort | uniq -i -c | grep " 2 " | cut -c9- >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred-and-feelnc.noncoding\
    grep "isBest\|^1" /output/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.best; head -n1 /output/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered; head -n1 /output/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best; head -n1 /output/hydra19H1_lncRNA-potential_FEELnc_DATE/stem/codport_stem/stem.lncRNAclasses.txt >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best.intergenic; awk 'NR==FNR{a[tolower($1)]; next} tolower($2) in a' /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.ezlncpred-and-feelnc.noncoding /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.best >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best; cat /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best | awk 'NR==1{print}; NR>1{if($7=="intergenic"){print}}' >> /output/hydra20H1_predicted-lncRNAs_Bash_${thislogdate}/transcripts/transcripts.feelncmain.filtered-and-best.intergenic\
    "

###############################################################
#                                                             #
#          STEP 21 [HYBRID]: RETRIEVE lncRNA METRICS          #
#                                                             #
###############################################################

## STEP 21H1 [resources: -m 1 -c 3 -w "01:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
docker run --rm \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    cat /output/hydra20H1_predicted-lncRNAs_Bash_DATE/transcripts_spliced/transcripts_spliced.ezlncpred-and-feelnc.noncoding | while read transcript; do grep "gene_id \"${transcript}\";" /output/hydra19H1_lncRNA-potential_FEELnc_${thislogdate}/transcripts_spliced/transcripts_spliced.gtf >> /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gtf; done\
    perl tools/gtf2gff3.pl /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gtf >> /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gff\
    /tools/genestats.sh /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.gff >> /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.genestats\
    cat /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.genestats | while read line; do echo "${line}" | awk -v OFS="\t" '{ exon_sum += $3; exon_len += $4; intron_sum += $5; intron_len += $6 } END { print exon_sum, exon_len / exon_sum, intron_sum, intron_len / (intron_sum+1) }' >> /output/hydra21H1_metrics-lncRNAs_GeneStats_${thislogdate}/transcripts_spliced/transcripts_spliced.metrics ; done\
    "

## STEP 22H1 [resources: -m 2 -c 20 -w "60:00:00"]
# WARNING 1: Reproduce for unspliced
# WARNING 2: Same step can be done for Trinity assembly. Simply replace input info and stem of outputs
# WARNING 3: User must provide:
reference_fasta_dir="/path/from/working/dir/to/ref-transcriptome_dir"
transgtf_dir="/path/from/working/dir/to/ref-transcriptome_dir"
lncrnadb_dir="/path/from/working/dir/to//data_lncrnadb_dir"
lncrnadb_fasta_dir="/path/from/working/dir/to/lncRNAdb.fasta_dir" #bedtools getfasta -name -fi genome.fa -bed /path/from/working/dir/to/lncRNAdb.bed > lncRNAdb.fasta

docker run --rm \
  -v "${reference_fasta_dir}":/data_reference_fasta_dir \
  -v "${transgtf_dir}":/data_transgtf_dir \
  -v "${lncrnadb_fasta_dir}":/data_lncrnadb_fasta_dir \
  -v "${project_output_dir}":/output \
  --entrypoint /bin/bash \
  hydra:1.0.0 \
  -c "\
    thislogdate=$(date +'%d%m%Y%H%M%S%Z') && \
    cut -f2 /output/hydra20H1_predicted-lncRNAs_Bash_DATE/transcripts_spliced/transcripts.spliced.feelncmain.filtered-and-best | sort | uniq >> /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.ids ; seqtk subseq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.singleline.fasta /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.ids > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.fasta\
    pblat -threads=20 -minIdentity=75 /data_reference_fasta_dir/ref-transcriptome.fasta /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.fasta /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.psl\
    pblat -threads=20 -minIdentity=75 /data_lncrnadb_fasta_dir/lncRNAdb.fasta /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.fasta /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.psl\
    awk -F'\t' '{if($1>=0.85*$11 && $1>=0.85*$15){print}}' /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.psl | sed -e '1,2d' > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.psl\
    cut -f14 /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.psl | while read t; do grep -P "\ttranscript\t.*${t}" /data_transgtf_dir/ref-transcriptome.gtf ; done | sort | uniq > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.gtf\
    cat /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.psl | while read l; do t=`echo ${l} | cut -d" " -f14`; type=`grep "${t}" /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.gtf | cut -d'"' -f20` ; echo -e "${l}\t${type}"; done > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.types\
    awk -F'\t' '{if($1>=0.85*$11 && $1>=0.85*$15){print}}' /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.psl | sed -e '1,2d' > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.85c.psl\
    cat /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.lncrnadb.85c.psl | while read l; do echo -e "${l}\tlncrnadb"; done >> /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.types\
    cat /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.types | while read l; do type=`echo ${l} | cut -d" " -f22`; id=`echo ${l} | cut -d" " -f10`; if [[ "$type" == "lincRNA" || "$type" == "lncrnadb" ]]; then echo -e "${id}" >> /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.lnc; else echo -e "${id}" >> /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc.tmp; fi; done ; sort /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc.tmp | uniq > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc; rm -f /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc.tmp\
    sort /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.ids /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.notlnc | uniq -c | grep " 1 " | tr -s ' ' | cut -d' ' -f3 > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc\
    cat /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc | while read t; do grep -P "\ttranscript\t.*${t}" /data_transgtf_dir/ref-transcriptome.gtf ; done > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc.gtf\
    seqtk subseq ${project_output_dir}/hydra11H1_assembly_rnaSPAdes_DATE/transcripts.singleline.fasta /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc > /output/hydra21H2_lncRNA-annotation_PBLAT_${thislogdate}/transcripts_spliced/transcripts_spliced.75i.ref-stem.85c.final.lnc.fasta\
    "

###############################################################
#                                                             #
#         STEP 22 [HYBRID]: SUMMARY lncRNA DISCOVERY          #
#                                                             #
###############################################################

## STEP 22H1 [resources: -m 1 -c 1 -w "24:00:00"]
# This step produces a series of TSV files from out and error files printed from stdout on a PBS. Please adapt where necessary to suit how you've run the steps so far.
# Please go through each  series of commands on write-to-pbs/ipda_HYDRA_step22H1-to-pbs.sh 