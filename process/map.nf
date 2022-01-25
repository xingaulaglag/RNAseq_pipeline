/*
*Scripts from Xingaulag
*Date: 10.01.2022 
*RNAseq pipeline from QC to mapping
*/

// nextflow.enable.dsl=2

//     params.reads="ercc_samples/samples/*{1,2}.fastq"
//     params.length=20
//     params.quality=20
//     params.genome="ercc_samples/ref/chr22_with_ERCC92.fa"
//     params.out_Dir="ercc_samples/output"
//     params.qc="$params.out_Dir/QC"
//     params.trim="$params.out_Dir/trim"
//     params.map="$params.out_Dir/map"
//     params.genomeDir='/mnt/d/RNAseq/ercc_samples/index'
        

// log.info """\
//          R N A S E Q - N F   P I P E L I N E    
//          ===================================
//          reads        : ${params.reads}
//          length       : ${params.length}
//          quality      : ${params.quality}
//          Ref genome   : ${params.genome}
//          Index genome : ${params.genomeDir}
//          Output storage: ${params.out_Dir}
//          """
//          .stripIndent()

//  Channel
//     .fromFilePairs(params.reads, checkIfExists: true)
//     .set{read_pairs_ch}
    
process qc {
    tag "QC from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/QC", mode:'copy'
    
    input:
    tuple val(pair_id), path(reads)

    output:
    path ('*_fastqc.{zip,html}'), emit: report

    script:
    """
    fastqc ${reads} --threads $task.cpus -q
    """
}
process trimming {
    tag "Trimming from $pair_id"
    publishDir "${baseDir}/ercc_samples/output/trim", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*fq"), emit: trim_out
    path '*trimming_report.txt', emit: report
    path '*_fastqc.{zip,html}', emit: report_trim_qc

    script:
    """
    trim_galore --paired --fastqc --retain_unpaired --length ${params.length} \\
    -j $task.cpus ${reads[0]} ${reads[1]}
    
    """
}

genome=Channel.fromPath(params.genome)
index=Channel.fromPath(params.genomeDir)


process mapping {

    tag "Mapping from $pair_id"
    publishDir "$params.map", mode: 'copy'
    
    input:
    tuple val(pair_id), path(trim)
    path genomeDir
    
    output:
    tuple val(pair_id), path("${pair_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    path("${pair_id}_Aligned.sortedByCoord.out.bam.bai"), emit:bai
    path("*.out"), emit: report


    script:
    """
    STAR --readFilesIn ${params.trim}/${pair_id}1_val_1.fq ${params.trim}/${pair_id}2_val_2.fq \\
    --genomeDir ${params.genomeDir} \\
    --outFileNamePrefix ${pair_id}_ \\
    --runThreadN ${task.cpus} \\
    --outSAMtype BAM SortedByCoordinate \\
    --peOverlapNbasesMin 10 \\
    --alignIntronMax 1 \\
    --peOverlapMMp 0.01

    samtools index ${pair_id}_Aligned.sortedByCoord.out.bam
    samtools index -c ${pair_id}_Aligned.sortedByCoord.out.bam
    """
}


workflow {

    qc(read_pairs_ch)

    trimming(read_pairs_ch)

    mapping(trimming.out.trim_out, params.genomeDir)
}