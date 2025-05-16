// modules/step01S1.nf
process step01S1 {
    tag "$f1"

    input:
    path f1
    path f2
    val logdate
    val taskName

    output:
    path "hydra01S1_qc-raw_FastQC_${logdate}_${taskName}/"

    container 'hydra/containers/hydra_nonconflicting:latest'

    cpus params.cpus
    memory params.memory
    time params.time

    script:
    def logdate = new Date().format("ddMMyyyyHHmmssZ")
    def outdir = "hydra01S1_qc-raw_FastQC_${logdate}_${task.runName}"
    """
    mkdir -p ${outdir}
    fastqc -t ${task.cpus} --outdir ${outdir} ${f1} ${f2} --memory ${task.memory.toMega()}M
    
    &

    if [[ $flag_merge = "yes" ]]; then zcat -c ${path_file}*_*1.f* | gzip -c >> ${out_path_step02S1_Rcorrector}/${file}_R1.fastq.gz ; fi
    if [[ $flag_merge = "yes" ]]; then zcat -c ${path_file}*_*2.f* | gzip -c >> ${out_path_step02S1_Rcorrector}/${file}_R2.fastq.gz ; fi
    if [[ $flag_merge = "yes" ]]; then f1="${out_path_step02S1_Rcorrector}/${file}_R1.fastq.gz"; f2="${out_path_step02S1_Rcorrector}/${file}_R2.fastq.gz"; elif [[ $flag_merge = "no" ]]; then f1=`grep "${path_file}" ${input} | cut -f2`; f2=`grep "${path_file}" ${input} | cut -f3`; fi; perl ${module_rcorrector} -t ${ncpus} -od ${out_path_step02S1_Rcorrector} -1 ${f1} -2 ${f2}
    f1="${out_path_step02S1_Rcorrector}/${file}_R1.cor.fq.gz"; f2="${out_path_step02S1_Rcorrector}/${file}_R2.cor.fq.gz"; bash ${module_reformat} in=${f1} in2=${f2}

    """
}
