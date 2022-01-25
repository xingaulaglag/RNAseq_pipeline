/*
*Scripts form Xingaulag
*Date: 17.01.2022 
*RNAseq pipeline for Qualimap only
*/
// nextflow.enable.dsl=2

// // baseDir ="/home/xingau/workflow"
// params.baseDir='/mnt/d/RNAseq'
// params.output="/home/xingau/workflow/ercc_samples/output/"
// params.bed12="/home/xingau/workflow/ercc_samples/ref/bed12/chr22_with_ERCC92.bed12"
// params.csvDir = '/home/xingau/workflow/metadata/metadata.csv' 
// params.gtf="/home/xingau/workflow/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"

// meta = Channel.from(file(params.csvDir))
//               .splitCsv(header:true)
//               .map{ row-> tuple("$row.pair_id"), file("$row.bam_path"), file("$row.bai_path")}
//               .set{sample_ch}

// log.info """\
//          R N A S E Q - N F   P I P E L I N E    
//          ===================================
//          Genome :           ${params.bed12}
//          Metadata:          ${params.csvDir}
//          Output storage:    ${params.output}
//          """
//          .stripIndent()

process qualimap {
    tag "Qualimap from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/qualimap" , mode:'copy'
    
    input:
    tuple val(pair_id), path(bam)
    path(params.gtf)

    output:
    path("${pair_id}_Aligned.sortedByCoord.out_rnaseq_qc/${pair_id}.pdf") , emit: report
    // file "v_qualimap.txt"

    script:
    """
    qualimap rnaseq -bam ${bam} -gtf ${params.gtf} -outfile ${pair_id}.pdf -p strand-specific-forward
    """
}
workflow {
 
   qualimap(mapping.out.bam, params.gtf)
}