/*
*Scripts form Xingaulag
*Date: 17.01.2022 
*RNAseq pipeline for RSeqc only
*/
// nextflow.enable.dsl=2
// params.baseDir='/mnt/d/RNAseq'
// params.output="${baseDir}/ercc_samples/output/"
// params.bed12="${baseDir}/ercc_samples/ref/bed12/chr22_with_ERCC92.bed12"
// // params.csvDir ='${baseDir}/metadata/ercc_metadata.csv' 
// params.gtf="${baseDir}/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"

// meta = Channel.from(file(params.csvDir))
//               .splitCsv(header:true)
//               .map{ row-> tuple("$row.pair_id"), file("$row.bam_path"), file("$row.bai_path") }
//               .set{sample_ch}

// log.info """\
//          R N A S E Q - N F   P I P E L I N E    
//          ===================================
//          Genome :           ${params.bed12}
//          Metadata:          ${params.csvDir}
//          Output storage:    ${params.output}
//          """
//          .stripIndent()

process rseqc {
    tag "RNA seq QC from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/rseqc" , mode:'copy'
    
    input:
    tuple val(pair_id), path(bam)
    path(bai)
    path(params.bed12)

    output:
    path "${pair_id}*" , emit: report

    script:
    """
    read_distribution.py -i ${bam} -r ${params.bed12} > "${pair_id}.rseqc.read_distribution.txt"
    infer_experiment.py -i ${bam} -r ${params.bed12}
    inner_distance.py -i ${bam} -o ${pair_id}.rseqc -r ${params.bed12}
    read_duplication.py -i ${bam} -o ${pair_id}.rseqc.read_duplication
    junction_annotation.py -i ${bam} -r ${params.bed12} -o ${pair_id}.rseqc.junction_annotation
    junction_saturation.py -i ${bam} -o ${pair_id}.rseqc -r ${params.bed12}
    bam_stat.py -i ${bam} > ${pair_id}.rseqc.bam_stat.txt

    """
}


workflow {
    rseqc(mapping.out.bam, mapping.out.bai, params.bed12)
}