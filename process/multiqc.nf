// // Generate report via MultiQC
// nextflow.enable.dsl=2
// params.publish_dir = "/home/xingau/workflow/ercc_samples/output/multiqc"

// params.skip_multiqc = false
// params.run_name = false
// params.ensembl_web = false
// params.ignore_R1 = false
// params.deseq2_fdr = 0.05
// params.gprofiler_fdr = 0.05
// params.cloudfront_link_duration = 60

process multiqc {

    publishDir "${baseDir}/ercc_samples/output/multiqc", mode: 'copy'

    input:
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path (report)
    path(report_trim_qc)

    output:
    file("*")

    script:
    
    """
    multiqc . -m fastqc -m Trim_Galore -m star -m rseqc -m preseq -m picard -m qualimap -m featureCounts \\
    -m custom_content -m plot_sample_distance -m plot_gene_heatmap -m DESeq2 -m gProfiler
    """
}

workflow {
    multiqc(qc.out.report.collect(), trimming.out.report, trimming.out.report_trim_qc, mapping.out.report.collect(), rseqc.out.report, mark_dup.out.report, preseq.out.report, dupradar.out.report, qualimap.out.report, featurecounts.out.report, deseq2.out.report, gprofiler.out.report)
}