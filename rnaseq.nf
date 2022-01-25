nextflow.enable.dsl=2

params.baseDir='/mnt/d/RNAseq'
    params.reads="${baseDir}/ercc_samples/samples/*{1,2}.fastq"
    params.length=20
    params.quality=20
    params.genome="${baseDir}/ercc_samples/ref/chr22_with_ERCC92.fa"
    params.out_Dir="${baseDir}/ercc_samples/output"
    params.qc="${baseDir}/ercc_samples/output/QC"
    params.trim="${baseDir}/ercc_samples/output/trim"
    params.map="${baseDir}/ercc_samples/output/map"
    params.genomeDir="${baseDir}/ercc_samples/index"
    params.bed12="${baseDir}/ercc_samples/ref/bed12/chr22_with_ERCC92.bed12"
    params.csvDir ="${baseDir}/metadata/ercc_metadata.csv"
    params.gtf="${baseDir}/ercc_samples/ref/gtf/chr22_with_ERCC92.gtf"
    params.bam_suffix = "_Build37-ErccTranscripts-chr22.read_Aligned.sortedByCoord.out.bam"
    params.design="${baseDir}/metadata/design.csv"
    params.compare="${baseDir}/metadata/comparison.csv"
    params.deseq2_fdr = 0.05
    params.deseq2_lfc = 0.585
    params.gprofiler_fdr = 0.05
    params.gprofiler_organism = 'hsapiens'

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set{read_pairs_ch}

// meta = Channel.from(file(params.csvDir))
//     .splitCsv(header:true)
//     .map{ row-> tuple("$row.pair_id"), file("$row.bam_path"), file("$row.bai_path") }
//     .set{sample_ch}

include { qc; trimming; mapping } from './process/map'
include { rseqc } from ('./process/rseqc')
include { qualimap } from ('./process/qualimap')
include { mark_dup; preseq } from ('./process/markdup_preseq')
include { dupradar } from ('./process/dupradar')
include { featurecounts; merge_counts } from ('./process/featurecounts_merge')
include { deseq2 } from ('./process/DESeq2.nf')
include { gprofiler } from ('./process/gprofiler')
include { multiqc } from ('./process/multiqc')

workflow {
    qc(read_pairs_ch)
    trimming(read_pairs_ch)
    mapping(trimming.out.trim_out, params.genomeDir)
    rseqc(mapping.out.bam, mapping.out.bai, params.bed12)
    qualimap(mapping.out.bam, params.gtf)
    mark_dup(mapping.out.bam, mapping.out.bai)
    preseq(mapping.out.bam, mapping.out.bai)
    dupradar(mapping.out.bam, mapping.out.bai, params.gtf)
    featurecounts(mapping.out.bam, mapping.out.bai, params.gtf)
    merge_counts(featurecounts.out.counts.collect())
    deseq2(merge_counts.out.merge, params.design, params.compare)
    gprofiler(deseq2.out.deseq2_result)
    
    multiqc(qc.out.report.collect(), \
    trimming.out.report.collect(), \
    trimming.out.report_trim_qc.collect(), \
    mapping.out.report.collect(), \
    rseqc.out.report.collect(), \
    mark_dup.out.report.collect(), \
    preseq.out.report.collect(), \
    dupradar.out.report.collect(),\
    qualimap.out.report.collect(), \
    featurecounts.out.report.collect(), \
    deseq2.out.report.collect(), \
    gprofiler.out.report.collect())
} 