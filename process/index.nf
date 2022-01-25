nextflow.enable.dsl=2

        // params.reads="$baseDir/ercc_samples/samples/*{1,2}.fastq.gz"
        // params.length=20
        // params.quality=20
        // params.genome="$baseDir/ercc_samples/ref/chr22_ERCC92_transcripts.fa"
        // params.out_Dir="$baseDir/ercc_samples/output"
        // params.qc="$params.out_Dir/QC"
        // params.trim="$params.out_Dir/trim"
        // params.map="$params.out_Dir/map"
        // params.index="$params.out_Dir/index"
        // params.genomeDir="$baseDir/ercc_samples/index"
        // params.bed12="$baseDir/ercc_samples/ref/bed12/chr22_rRNA.bed12"

genome=Channel.fromPath(params.genome)
index=Channel.fromPath(params.genomeDir)

process index {

    tag "Index ref Genome"
    publishDir "$params.index", mode:'copy'

    input:
    path(genome)

    output:

    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeFastaFiles  /home/xingau/workflow/ercc_samples/ref/chr22_with_ERCC92.fa --genomeDir /home/xingau/workflow/ercc_samples/index --genomeSAindexNbases 10
    """
}
workflow {

    index(genome)

}
