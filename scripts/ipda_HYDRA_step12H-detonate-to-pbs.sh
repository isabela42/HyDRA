#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Oct 04, 2023
Last modified on Oct 20, 2023
Version: ${version}

Description: Write and submit PBS jobs for Step 12H4 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step12H4-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/HYDRA_step11H*_de-novo-transcriptome-assembly_*_DATE/stem.fasta
                            of either transcripts.fasta or .Trinity.fasta files.

                            Col2:
                            path/from/working/dir/to/HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R1.cor.fq.gz
                            of all header-fixed files that were used to create the assembly,
                            seppareted by comma (same order as Col3)

                            Col3:
                            path/from/working/dir/to/HYDRA_step06S1_ribodepletion_BBDuk_DATE/stem_R2.cor.fq.gz
                            of all header-fixed files that were used to create the assembly,
                            seppareted by comma (same order as Col2)

                            Col4:
                            path/from/working/dir/to/HYDRA_step06L1_ribodepletion_BBDuk_DATE/stem.fasta
                            of all long-read derived fasta files that were then used to create the hybrid assembly,
                            seppareted by comma. Leave it empty for short-reads only assembly
                            

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

#   HYBRID DE NOVO TRANSCRIPTOME ASSEMBLY
#   ------------------------------------------------------------------------------------------------------------
#   12 [H] Raw assembly quality check (1BUSCO, 2rnaQUAST, 3HISAT2, 4DETONATE, 5TrinityStats) 

This tool could not be included in the final HYDRA pipeline version.

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
        ?) echo script usage: bash ipda_HYDRA_step12H4-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step12H4-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

reference_fasta=/working/lab_julietF/isabelaA/data/references/TranscriptomeGRCh38rel79.fasta
referencedir=/working/lab_julietF/isabelaA/data/references/
reference_index=/working/lab_julietF/isabelaA/data/references/TranscriptomeGRCh38rel79_bowtie2idx

#................................................
#  Required modules, softwares and libraries
#................................................

## Load tools from HPC
# For more info, see
# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>

# DETONATE
# <https://github.com/deweylab/detonate>
module_detonate=conda-envs/detonate

# UCSC Tools 20160223 ( BLAT 36x1)
module_ucsctools=ucsctools/20160223

# PBLAT
module_pblat=conda-envs/pblat-2.5.1

# GMAP 2023-07-20
module_gmap=gmap/20230720

# Bowtie2
module_bowtie2=bowtie2/2.2.9

# SamTools
module_samtools=samtools/1.9

# RSEM 1.3.1
module_rsem=RSEM/1.3.1

## Path to user-installed tools

# None required

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step12H4_DETONATE="HYDRA_step12H4_raw-assembly-quality-check_DETONATE_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step12H4_DETONATE}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step12H4-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step12H4_DETONATE}"
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
echo "## Executing bash ipda_HYDRA_step12H4-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step12H4_DETONATE}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#!/bin/sh" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Script:  ${pbs_stem}_${file}_${thislogdate}.pbs" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Version: v01" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Email:   ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#PBS -N ${pbs_stem}_${file}_${thislogdate}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#PBS -r n" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#PBS -m abe" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#PBS -M ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Set main working directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "## Change to main directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Load tools from HPC"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## For more info, see"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load DETONATE module at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_detonate}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load PBLAT at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_pblat}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load UCSC Tools 20160223 (BLAT 36x1) at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_ucsctools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Verify BLAT 36x1 version at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'blat | head -n1' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load GMAP 2023-07-20 module at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_gmap}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Verify GMAP 2023-07-20 version at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'gmap --version' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load Bowtie2 2.2.9 module at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_bowtie2}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Verify Bowtie2 version at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'bowtie2 --version' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load SamTools 1.9 module at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_samtools}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Verify SamTools version at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'samtools --version' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Load RSEM 1.3.1 module at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_rsem}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# Verify RSEM 1.3.1 version at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'rsem-calculate-expression --version' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Path to user-installed tools"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# None required"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write to retrieve working directory information
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Info on working directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Report file system disk space usage at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "df -H" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "echo" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Report files on working directory at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "ls -lh" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "echo" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run BLAT alignment at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; ref=`echo "$(basename "${reference_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; mkdir -p ${out_path_step12H4_DETONATE}/${file}/ ; echo "pblat -threads=${ncpus} -minIdentity=80 ${reference_fasta} ${path_file} ${out_path_step12H4_DETONATE}/${file}/${file}_to_${ref}.psl" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; ref=`echo "$(basename "${reference_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; echo "pblat -threads=${ncpus} -minIdentity=80 ${path_file} ${reference_fasta} ${out_path_step12H4_DETONATE}/${file}/${ref}_to_${file}.psl" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run GMAP alignment at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; ref=`echo "$(basename "${reference_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; longreads_fasta=`grep ${path_file} ${input} | cut -f4` ; [ -z "$longreads_fasta" ] && echo "File not provided" || echo "gmap -n 0 -t ${ncpus} -B 5 --dir=${referencedir} --gseg ${reference_fasta} ${longreads_fasta} | samtools view -@${ncpus} -bS - >> ${out_path_step12H4_DETONATE}/${file}/${file}_longreads_gmap_${ref}.bam" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Build reference bowtie2 index at at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; ref=`echo "$(basename "${reference_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; echo "bowtie2-build --threads ${ncpus} ${reference_fasta} ${reference_index}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

FAILED FROM HERE

cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run Bowtie2 alignment at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; f1=`grep ${path_file} ${input} | cut -f2` ; f2="grep ${path_file} ${input} | cut -f3"; echo "bowtie2 --fast --no-unal --threads ${ncpus} -x ${reference_index} -1 ${f1} -2 ${f2} | samtools view -@${ncpus} -bS - >> ${out_path_step12H4_DETONATE}/${file}/${file}_shortreads_bowtie2_${ref}.bam" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run Samtools merge at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "if [ -f ${out_path_step12H4_DETONATE}/${file}/${file}_longreads_gmap_${ref}.bam ]; then samtools merge -o ${out_path_step12H4_DETONATE}/${file}/${file}_merged_shortreads_bowtie2_and_longreads_gmap_${ref}.bam -@${ncpus} ${out_path_step12H4_DETONATE}/${file}/${file}_shortreads_bowtie2_${ref}.bam ${out_path_step12H4_DETONATE}/${file}/${file}_longreads_gmap_${ref}.bam ; bam=${out_path_step12H4_DETONATE}/${file}/${file}_longreads_gmap_${ref}.bam ; else echo 'No files to merge'; bam=${out_path_step12H4_DETONATE}/${file}/${file}_shortreads_bowtie2_${ref}.bam; fi" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Coordination sort transcriptome alignment at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "samtools sort -@${ncpus} ${bam} -o ${out_path_step12H4_DETONATE}/${file}/${file}_merged_shortreads_bowtie2_and_longreads_gmap_${ref}.coordsorted.bam" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run RSEM prepare expression at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; ref=`echo "$(basename "${reference_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; echo "rsem-prepare-reference --num-threads ${ncpus} ${reference_fasta} ${ref}_rsemindex-ref" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "rsem-prepare-reference --num-threads ${ncpus} ${path_file} ${file}_rsemindex-ref" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run RSEM calculate expression at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; ref=`echo "$(basename "${reference_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`; echo "rsem-calculate-expression --num-threads ${ncpus} --no-bam-output --alignments --paired-end ${bam} ${ref}_rsemindex-ref ${ref}_rsem-expr" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "rsem-calculate-expression --num-threads ${ncpus} --no-bam-output --alignments --paired-end ${bam} ${file}_rsemindex-ref ${file}_rsem-expr" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Move RSEM files to out folder at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "mv "${ref}_rsemindex-ref"* "${file}_rsemindex-ref"* "${ref}_rsem-expr"* "${file}_rsem-expr"* ${out_path_step12H4_DETONATE}/${file}/" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run DETONATE ref-eval at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; len=`awk '/^>/ {if (seq != "") {if (min_length == 0 || length(seq) < min_length) {min_length = length(seq)}; seq = ""} next} {seq = seq $0} END {if (min_length == 0 || length(seq) < min_length) {min_length = length(seq)}; print min_length}' ${path_file}` ; echo "ref-eval --scores=nucl,pair,contig,kmer,kc --weighted=both --A-seqs ${path_file} --B-seqs ${reference_fasta} --A-expr ${out_path_step12H4_DETONATE}/${file}/${file}_rsem-expr.isoforms.results --B-expr ${out_path_step12H4_DETONATE}/${file}/${ref}_rsem-expr.isoforms.results --A-to-B ${out_path_step12H4_DETONATE}/${file}/${file}_to_${ref}.psl --B-to-A ${out_path_step12H4_DETONATE}/${file}/${ref}_to_${file}.psl --num-reads 5000000 --readlen ${len} --kmerlen ${len} | tee ${out_path_step12H4_DETONATE}/${file}/detonate_${file}_scores.txt" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n31 ${pbs} | head -n1)" ; echo "#                       $(tail -n30 ${pbs} | head -n1)" ; echo "#                       $(tail -n27 ${pbs} | head -n1)" ; echo "#                       $(tail -n24 ${pbs} | head -n1)" ;  echo "#                       $(tail -n21 ${pbs} | head -n1)" ; echo "#                       $(tail -n18 ${pbs} | head -n1)" ; echo "#                       $(tail -n15 ${pbs} | head -n1)" ; echo "#                       $(tail -n12 ${pbs} | head -n1)" ; echo "#                       $(tail -n11 ${pbs} | head -n1)" ; echo "#                       $(tail -n8 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 12H4 jobs) at
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
sed -i 's,${module_detonate},'"${module_detonate}"',g' "$logfile"
sed -i 's,${module_pblat},'"${module_pblat}"',g' "$logfile"
sed -i 's,${module_ucsctools},'"${module_ucsctools}"',g' "$logfile"
sed -i 's,${module_rsem},'"${module_rsem}"',g' "$logfile"
sed -i 's,${reference_fasta},'"${reference_fasta}"',g' "$logfile"
sed -i 's,${referencedir},'"${referencedir}"',g' "$logfile"
sed -i 's,${reference_index},'"${reference_index}"',g' "$logfile"
sed -i 's,${out_path_step12H4_DETONATE},'"${out_path_step12H4_DETONATE}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v