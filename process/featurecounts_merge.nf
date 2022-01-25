nextflow.enable.dsl=2

// params.output="/home/xingau/workflow/ercc_samples/output/"
// params.csvDir = '/home/xingau/workflow/metadata/ercc_meta.csv' 
// params.gtf="/home/xingau/workflow/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"
// params.bam_suffix = "_Build37-ErccTranscripts-chr22.read_Aligned.sortedByCoord.out.bam"

// meta = Channel.from(file(params.csvDir))
//               .splitCsv(header:true)
//               .map{ row-> tuple("$row.pair_id"), file("$row.bam_path"), file("$row.bai_path")}
//               .set{sample_ch}

process featurecounts {
    tag "FeatureCounts from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/featurecounts" , mode:'copy'

    input:
    tuple val(pair_id), path(bam)
    path(bai)
    path params.gtf

    output:
    path("${pair_id}_gene.featureCounts.txt"), emit: counts
    path "${pair_id}_gene.featureCounts.txt.summary", emit:report
    // path "${pair_id}_transcript_id.featureCounts.txt"
    // path "${pair_id}_transcript_id.featureCounts.txt.summary"
    // path "v_featurecounts.txt"

    script:
    """
    featureCounts -a ${params.gtf} -F 'GTF' -g 'gene_id' -t 'exon' -o ${pair_id}_gene.featureCounts.txt --extraAttributes 'gene_name' -p -s 1 ${bam} -T 1
    """
}

process merge_counts {
    publishDir "${baseDir}/ercc_samples/output/merge_counts", mode: 'copy'

    input:
    path(counts)

    output:
    path 'merged_gene_counts.txt', emit:merge
    file 'gene_lengths.txt'

    script:

    gene_ids = "<(tail -n +2 ${counts[0]} | cut -f1,7 )"
    merge_counts = counts.collect{filename ->
    "<(tail -n +2 ${filename} | sed 's:${params.bam_suffix}::' | cut -f8)"}.join(" ")
    """
    paste $gene_ids $merge_counts > merged_gene_counts.txt
    tail -n +2 ${counts[0]} | cut -f1,6 > gene_lengths.txt
    """
}
workflow {
 
   featurecounts(mapping.out.bam, mapping.out.bai, params.gtf)
   merge_counts(featurecounts.out.counts.collect())
}