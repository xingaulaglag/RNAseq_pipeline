// Run pathway enrichment analysis with gProfiler
// nextflow.enable.dsl=2
// params.publish_dir = "/home/xingau/workflow/ercc_samples/output/gProfiler"
// params.gprofiler_organism = 'hsapiens'
// params.deseq2_fdr = 0.05
// params.gprofiler_fdr = 0.05
// params.deseq_results="/home/xingau/workflow/ercc_samples/output/deseq2/HBR_vs_UHR_DESeq_results.tsv"

process gprofiler {
    publishDir "${baseDir}/ercc_samples/output/gprofiler", mode: 'copy'

    input:
    path(deseq2_result)

    output:
    path "*_gProfiler_results.tsv", emit: report
    path "*_gProfiler_results.xlsx"
    path "v_gProfiler.txt"

    script:
    """
    gProfiler.py ${deseq2_result} -o $params.gprofiler_organism -q $params.deseq2_fdr -p $params.gprofiler_fdr
    pip freeze | grep gprofiler > v_gProfiler.txt
    """
}

workflow {

    gprofiller(deseq2.out.deseq2_result)
}