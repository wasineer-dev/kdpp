#!/usr/bin/env nextflow

// try 3.9075, 1.2, 0.1
params.threshold = 1.2
params.bootstrapp = "bootstrapp"
params.ndup = 1
params.samplerate = 0.9
params.sampletarget = 320000000

process sub_sample {
  cache false  
  maxForks 1
  publishDir "$params.bootstrapp/$sampleId"

  input:
  tuple val(rnd_seed), val(sampleId), path(bbduk_files)

  output:
  tuple val({rnd_seed}), val({sampleId}), path("*sub*.fq")
  
  """
    seqtk sample -s${rnd_seed} ${bbduk_files[0]} 0.5 > ${sampleId}_${rnd_seed}_sub1.fq
    seqtk sample -s${rnd_seed} ${bbduk_files[1]} 0.5 > ${sampleId}_${rnd_seed}_sub2.fq
  """
}

// Not pararellize
process build_vlmc {
    cache false
    maxForks 4
    publishDir "$params.bootstrapp/$sampleId"
    
    input:
    tuple val(rnd_seed), val(sampleId), path(bbduk_files)
    
    output:
    tuple val({rnd_seed}), path("*.bintree")

    """
        reformat.sh in1=${bbduk_files[0]} in2=${bbduk_files[1]} \
            out1=${sampleId}_${rnd_seed}_sub1.fq.gz out2=${sampleId}_${rnd_seed}_sub2.fq.gz \
            samplerate=$params.samplerate \
            samplebasestarget=$params.sampletarget

        cat ${sampleId}_${rnd_seed}_sub1.fq.gz ${sampleId}_${rnd_seed}_sub2.fq.gz > ${sampleId}_${rnd_seed}_merged.fq.gz
        dvstar \
        -m build \
        -p ${sampleId}_${rnd_seed}_merged.fq.gz \
        --out-path ${sampleId}_${rnd_seed}.bintree \
        --threshold $params.threshold --max-depth 4 --min-count 100 \
        --adjust-for-sequencing-errors \
        --sequencing-depth 20.0 \
        --sequencing-error-rate 0.001 \
        --pseudo-count-amount 0.01 \
        --temp-path tmp
        rm -rf ${sampleId}_${rnd_seed}_sub1.fq.gz ${sampleId}_${rnd_seed}_sub2.fq.gz ${sampleId}_${rnd_seed}_merged.fq.gz
    """
}

process VLMCDistance {
    cache false
    publishDir 'dvstar_distance'

    input:
    tuple val(rndId), path(source), path(target, stageAs: "?/*")

    output:
    stdout

    shell:
    """
        #!/bin/bash
        
        DISTANCE=`dvstar \
                        --mode dissimilarity --dissimilarity dvstar \
                        --in-path !{source} \
                        --to-path !{target} \
                        | tail -n 1` 
        
        echo "!{rndId}\t!{source}\t!{target}\t\$DISTANCE" 2>/dev/null
    """
}

process nj_tree {
    input:
    path filename

    output:
    stdout

    """
        python /Users/wasinee/nextflow/create_njtree.py ${filename}
    """
}

rnd_ch = Channel.of(1000..1000000)
          .randomSample(params.ndup)
          
reads_ch = Channel.fromFilePairs("/Users/wasinee/nextflow/postTrim/S*/S*{R1,R2}_*.bbduk.fastq")
        .ifEmpty { error "Cannot find any reads matching" }
        // Sort the channel elements based on the first object of each tuple,
        // that is, the sample name, and convert to a channel with a single
        // element which is a list of tuples
        .toSortedList( { a, b -> a[0] <=> b[0] } ) // <=> is an operator for comparison
        // flatten the single-element channel to a channel with as many elements
        // as there are samples, which is the original structure provided by
        // fromFilePairs
        .flatMap()
 
workflow {

    source_ch = rnd_ch.combine( reads_ch )
                | build_vlmc 

    distance_ch = source_ch.combine( source_ch, by: 0)
        | VLMCDistance 
        | collectFile(storeDir:'./vlmc_distances', name: "${params.output}")
        | view
        //| nj_tree
        //| collectFile(name: './phylip/outtrees', newLine: false)
        //| view
}