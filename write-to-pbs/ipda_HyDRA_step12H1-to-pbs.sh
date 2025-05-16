#!/bin/bash

version="1.0.1"
usage(){
echo "
Written by Isabela Almeida
Created on Sep 26, 2023
Last modified on May 15, 2025
Version: ${version}

Description: Write and submit PBS jobs for Step 12H1 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HyDRA_step12H1-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline development: -m 2 -c 5 -w "50:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/hydra11H*_assembly_*_DATE/stem.fasta
                            of either transcripts.fasta or .Trinity.fasta files.

                            Col2:
                            path/from/working/dir/to/BuscoLineages/busco_odb9

                            Col3:
                            stem-of-analysis
                            e.g.: eukaryota-odb9-transcripts

                            Col4:
                            species
                            e.g.: human 
                            Note: See BUSCO's Augustus documentation for available options.

                            It does not matter if same stem.fasta 
                            appears more than once on this input file.
                            However, each cell in col3 must be unique.
                   
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
#-->12 [H] Raw assembly quality check (1BUSCO, 2TrinityStats, 3TransRate, 4TransRate)
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
        ?) echo script usage: bash ipda_HyDRA_step12H1-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_hydra12H1-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# Python3
module_python3=python/3.6.1

# BUSCO (Benchmarking Universal Single-Copy Orthologs)
module_busco=busco/20161119

# BLAST
module_oldblast=ncbi-blast/2.2.31

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step12H1_BUSCO="hydra12H1_qc_BUSCO_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step12H1_BUSCO}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HyDRA_step12H1-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input files:                 ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${out_path_step12H1_BUSCO}"
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
echo "## Executing bash ipda_HyDRA_step12H1-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step12H1_BUSCO}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#!/bin/sh" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "##########################################################################" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Script:  ${pbs_stem}_${stem}_${thislogdate}.pbs" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Version: v01" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Email:   ${email}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "##########################################################################" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write PBS directives
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#PBS -N ${pbs_stem}_${stem}_${thislogdate}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#PBS -r n" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#PBS -m abe" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#PBS -M ${email}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write directory setting
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Set main working directory" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "## Change to main directory" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write load modules
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "module load ${module_python3}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "module load ${module_busco}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "module load ${module_oldblast}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#  Run step" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo 'echo "## Run BUSCO at" ; date ; echo' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; busco_odb9=`grep ${stem} ${input} | cut -f2 | sort | uniq`;  species=`grep ${stem} ${input} | cut -f4 | sort | uniq`; echo "BUSCO.py -i ${path_file} -l ${busco_odb9} -o ${stem} -m transcriptome -c ${ncpus} -t ${out_path_step12H1_BUSCO} -sp ${species}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo 'echo "## Move BUSCO files to output folder at" ; date ; echo' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f3 ${input} | sort | uniq | while read stem; do path_file=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; echo "mv run_${stem} ${out_path_step12H1_BUSCO}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 12H1 jobs) at
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
sed -i 's,${module_python3},'"${module_python3}"',g' "$logfile"
sed -i 's,${module_busco},'"${module_busco}"',g' "$logfile"
sed -i 's,${module_oldblast},'"${module_oldblast}"',g' "$logfile"
sed -i 's,${species},'"${species}"',g' "$logfile"
sed -i 's,${busco_odb9},'"${busco_odb9}"',g' "$logfile"
sed -i 's,${stem_db},'"${stem_db}"',g' "$logfile"
sed -i 's,${out_path_step12H1_BUSCO},'"${out_path_step12H1_BUSCO}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v
