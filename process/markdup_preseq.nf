/*
*Scripts form Xingaulag
*Date: 17.01.2022 
*RNAseq pipeline for Mark Duplicates and Index Mark_dup_bam
*/
nextflow.enable.dsl=2

// params.output="/home/xingau/workflow/ercc_samples/output/"
// params.csvDir = '/home/xingau/workflow/metadata/metadata.csv' 
// params.gtf="/home/xingau/workflow/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"

// meta = Channel.from(file(params.csvDir))
//               .splitCsv(header:true)
//               .map{ row-> tuple("$row.pair_id"), file("$row.bam_path"), file("$row.bai_path")}
//               .set{sample_ch}

// log.info """\
//          R N A S E Q - N F   P I P E L I N E    
//          ===================================
//          Metadata:          ${params.csvDir}
//          Output storage:    ${params.output}
//          """
//          .stripIndent()
process mark_dup {
    tag "Qualimap from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/markdup" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    

    output:
    tuple val(pair_id), path("${pair_id}.markDups.bam"), path("${pair_id}.markDups.bam.{bai,csi}")
    path "${pair_id}*", emit: report

    script:
    """
    picard MarkDuplicates \\
            INPUT=${bam} \\
            OUTPUT=${pair_id}.markDups.bam \\
            METRICS_FILE=${pair_id}.markDups_metrics.txt \\
            REMOVE_DUPLICATES=false \\
            ASSUME_SORTED=true \\
            PROGRAM_RECORD_ID='null' \\
            VALIDATION_STRINGENCY=LENIENT
        samtools index -c ${pair_id}.markDups.bam
        samtools index ${pair_id}.markDups.bam 
    """
}
process preseq {
    tag "Qualimap from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/preseq" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    
    output:
    path "${pair_id}*", emit: report

    script:
    """
    preseq lc_extrap -v -B ${bam} -o ${pair_id}.ccurve.txt
    """
}
workflow {
 
   mark_dup(mapping.out.bam, mapping.out.bai)
   preseq(mapping.out.bam, mapping.out.bai)
}