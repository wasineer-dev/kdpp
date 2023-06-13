#!/usr/bin/env nextflow

// try 3.9075, 1.2, 0.1
params.threshold = 1.2

process BuildVLMC {
    publishDir 'dvstar_output/', pattern: "*.bintree"
    
    input:
    tuple val(sampleId), path(bbduk_files)
    
    output:
    tuple val({sampleId}), path("*.bintree")

    """
        reformat.sh in1=${bbduk_files[0]} in2=${bbduk_files[1]} \
            out1=${sampleId}_sub1.fq.gz out2=${sampleId}_sub2.fq.gz samplerate=0.8

        cat ${sampleId}_sub1.fq.gz ${sampleId}_sub2.fq.gz > ${sampleId}_merged.fq.gz
        dvstar \
        -m build \
        -p ${sampleId}_merged.fq.gz \
        --out-path ${sampleId}.bintree \
        --threshold $params.threshold --max-depth 4 --min-count 100 \
        --adjust-for-sequencing-errors \
        --sequencing-depth 20.0 \
        --sequencing-error-rate 0.001 \
        --pseudo-count-amount 0.01 \
        --temp-path tmp
    """
}

process vlmc_scores {
    //publishDir 'dvstar_output/scores', pattern: "*_scores.txt"

    input:
    tuple val(id), path(vlmc_file), val(sampleId), path(bbduk_files)
    
    output:
    stdout

    shell:
    """
        #!/bin/bash
        SCORE=`dvstar -m score \
                --fasta-path !{bbduk_files[0]} \
                --in-path !{vlmc_file} --max-depth 4 | tail -n 1`
        
        echo "!{id}\t!{sampleId}\t\$SCORE"
    """
}

bbduk_ch = Channel.fromFilePairs("/Users/wasinee/nextflow/postTrim/S*/S*{R1,R2}_*.bbduk.fastq")
        .ifEmpty { error "Cannot find any reads matching" }
        

workflow {
    vlmc_ch = Channel.fromFilePairs("/Users/wasinee/nextflow/postTrim/S*/S*{R1,R2}_*.bbduk.fastq")
                .ifEmpty { error "Cannot find any reads matching" }
                // Sort the channel elements based on the first object of each tuple,
                // that is, the sample name, and convert to a channel with a single
                // element which is a list of tuples
                .toSortedList( { a, b -> a[0] <=> b[0] } ) // <=> is an operator for comparison
                // flatten the single-element channel to a channel with as many elements
                // as there are samples, which is the original structure provided by
                // fromFilePairs
                .flatMap()
                .view()
                | BuildVLMC
                | combine(bbduk_ch)
                | view
                | vlmc_scores
                | collectFile(name: "./dvstar_output/scores.tsv", newLine: true)           
}