// Run study level analysis/plots and DEG analysis with DESeq2
nextflow.enable.dsl=2

// params.output="/home/xingau/workflow/ercc_samples/output/"
// params.merge_count="/home/xingau/workflow/ercc_samples/output/merge_counts/merged_gene_counts_tab.txt"
// params.design="/home/xingau/workflow/metadata/design.csv"
// params.compare="/home/xingau/workflow/metadata/comparison.csv"
// params.deseq2_fdr = 0.05
// params.deseq2_lfc = 0.585

// design = Channel.fromPath(params.design)
// merge_count = Channel.fromPath(params.merge_count)
// compare = Channel.fromPath(params.compare)

process deseq2 {
    tag "DESeq2"
    publishDir "${baseDir}/ercc_samples/output/deseq2", mode: 'copy'

    input:
    path(merge)
    path(params.design)
    path(params.compare)

    output:
    path "*.{xlsx,jpg}"
    path "*_DESeq_results.tsv", emit: deseq2_result
    path "*{heatmap,plot,matrix}.tsv", emit: report
    path "v_DESeq2.txt"

    script:
    """
    DESeq2.r ${merge} ${params.design} $params.deseq2_fdr $params.deseq2_lfc ${params.compare}
    Rscript -e "write(x=as.character(packageVersion('DESeq2')), file='v_DESeq2.txt')"
    """
}

workflow {
    deseq2(merge_counts.out.merge, params.design, params.compare)
}
