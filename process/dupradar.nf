/*
*Scripts form Xingaulag
*Date: 17.01.2022 
*RNAseq pipeline for Mark Duplicates and Index Mark_dup_bam
*/
// nextflow.enable.dsl=2

// params.baseDir='/mnt/d/RNAseq'

// params.output="${baseDir}/ercc_samples/output/"
// params.csvDir = '${baseDir}/metadata/ercc_metadata.csv' 
// params.gtf="${baseDir}/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"

// meta = Channel.from(file(params.csvDir))
//               .splitCsv(header:true)
//               .map{ row-> tuple("$row.pair_id"), file("$row.bam_path"), file("$row.bai_path")}
//               .set{sample_ch}
              
process dupradar {
    tag "Dupradar from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/dupradar" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    path params.gtf

    output:
    path "*.pdf"
    path "*_dupMatrix.txt"
    path "*_intercept_slope.txt"
    path "*_mqc.txt", emit: report
    path "v_dupRadar.txt"

    script:
    """
    dupRadar.r ${bam} ${params.gtf} 0 paired 1
    Rscript -e "write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
    """
}
workflow {
   dupradar(mapping.out.bam, mapping.out.bai, params.gtf)
}