#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Created on Oct 03, 2023
Last modified on Feb 06, 2024
Version: ${version}

Description: Write and submit PBS jobs for Step 12H2 (hybrid) of the
HYDRA pipeline (HYbrid De novo RNA Assembly pipeline).

Usage: bash ipda_HYDRA_step12H2-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:                            

                            Col1:
                            path/from/working/dir/to/HYDRA_step11H*_de-novo-transcriptome-assembly_*_DATE/stem.fasta
                            of either transcripts.fasta or .Trinity.fasta files.

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
        ?) echo script usage: bash ipda_HYDRA_step12H2-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_HYDRA_step12H2-to-pbs_${thislogdate}.txt

#................................................
#  Additional information
#................................................

reference_fasta=/working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_T2T009914755.4_transcriptome.fa
reference_gtf=/working/lab_julietF/isabelaA/data/references/T2TCHM13v2/ensembl_GCA_009914755.4-2022_07-transcripts.gtf

#................................................
#  Required modules, softwares and libraries
#................................................

## Load tools from HPC
# For more info, see
# <https://genomeinfo.qimrberghofer.edu.au/wiki/HPC/Avalon#Loading_Software_.28modules.29>

# Python3
module_python3=python/3.6.1

# GMAP 2023-07-20
module_gmap=gmap/20230720

# BLAST
module_blast=ncbi-blast/2.10.0+

# rnaQUAST 2.2.3
# <https://cab.spbu.ru/files/rnaquast/release1.5.2/manual.html>
module_rnaQUAST=conda-envs/rnaquast

## Path to user-installed tools

# None required

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
out_path_step12H2_rnaQUAST="HYDRA_step12H2_raw-assembly-quality-check_rnaQUAST_${thislogdate}"

## Create output directories
mkdir -p ${out_path_step12H2_rnaQUAST}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_HYDRA_step12H2-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step12H2_rnaQUAST}"
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
echo "## Executing bash ipda_HYDRA_step12H2-to-pbs.sh"
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
echo "## Output files saved to:       ${out_path_step12H2_rnaQUAST}"
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
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_python3}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_gmap}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_blast}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "module load ${module_rnaQUAST}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Path to user-installed tools"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "# None required"' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Run rnaQUAST at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; mkdir -p ${out_path_step12H2_rnaQUAST}/${file}/ ; echo "rnaQUAST.py -r ${reference_fasta} --gtf ${reference_gtf} -c ${path_file} -t ${ncpus} -o ${out_path_step12H2_rnaQUAST}/${file}/" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo 'echo "## Deactivate conda environment at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo "$(basename "${path_file%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'` ; echo "conda deactivate" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n4 ${pbs} | head -n1)" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including HYDRA step 12H2 jobs) at
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
sed -i 's,${module_gmap},'"${module_gmap}"',g' "$logfile"
sed -i 's,${module_blast},'"${module_blast}"',g' "$logfile"
sed -i 's,${module_rnaQUAST},'"${module_rnaQUAST}"',g' "$logfile"
sed -i 's,${reference_fasta},'"${reference_fasta}"',g' "$logfile"
sed -i 's,${reference_gtf},'"${reference_gtf}"',g' "$logfile"
sed -i 's,${out_path_step12H2_rnaQUAST},'"${out_path_step12H2_rnaQUAST}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v
