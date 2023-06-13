#!/usr/bin/env nextflow

// try 3.9075, 1.2, 0.1
params.threshold = 1.2
params.ndup = 2

process BuildVLMC {
    maxForks 4
    publishDir 'dvstar_output/', pattern: "*.bintree"
    
    input:
    tuple val(id), val(sampleId), path(bbduk_files)
    
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
        rm -rf ${sampleId}_sub1.fq.gz ${sampleId}_sub2.fq.gz ${sampleId}_merged.fq.gz
    """
}

process vlmc_scores {
    //publishDir 'dvstar_output/scores', pattern: "*_scores.txt"

    input:
    tuple val(id), path(bbduk_files), val(sampleId), path(vlmc_file)
    
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

rnd_ch = Channel.of(1000..1000000)
          .randomSample(params.ndup)

ids_ch = Channel.fromFilePairs("/Users/wasinee/nextflow/postTrim/S*/S*{R1,R2}_*.bbduk.fastq")
        .ifEmpty { error "Cannot find any reads matching" }
        .map{ id, file -> tuple(id, "${id}.bintree") }
        //.toSortedList( { a, b -> a[0] <=> b[0] } ) 
        
        
bbduk_ch = Channel.fromFilePairs("/Users/wasinee/nextflow/postTrim/S*/S*{R1,R2}_*.bbduk.fastq")
        .ifEmpty { error "Cannot find any reads matching" }
        .toSortedList( { a, b -> a[0] <=> b[0] } ) 
        .flatMap()

workflow {
    //bbduk_ch.combine(ids_ch)
    //    .view()

    vlmc_ch = rnd_ch.combine( bbduk_ch )
                | BuildVLMC
                
    bbduk_ch.combine(vlmc_ch)
            | vlmc_scores
            | collectFile(name: "./dvstar_output/scores.tsv", newLine: false)           
}