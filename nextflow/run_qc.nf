#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.OUTPUT = "multiqc_output"

process FastQC {
    input:
    tuple val(name), file(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
        fastqc $reads
    """
}

process Multiqc {
    // where to store the results and in which way
    publishDir(params.OUTPUT, mode: 'copy')

    input:
    path (inputfiles)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {
    Channel.fromFilePairs("/Users/wasinee/nextflow/postTrim/S*/S*_{R1,R2}_*.bbduk.fastq")
        .ifEmpty { error "Cannot find any reads matching" }
        // Sort the channel elements based on the first object of each tuple,
        // that is, the sample name, and convert to a channel with a single
        // element which is a list of tuples
        .toSortedList( { a, b -> a[0] <=> b[0] } ) // <=> is an operator for comparison
        // flatten the single-element channel to a channel with as many elements
        // as there are samples, which is the original structure provided by
        // fromFilePairs
        .flatMap()
        // View the channel elements by printing it to the screen
        .view()
        | FastQC | collect | Multiqc
}
