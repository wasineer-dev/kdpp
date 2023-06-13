#!/usr/bin/env nextflow

// try 3.9075, 1.2, 0.1
params.threshold = 1.2
params.bootstrapp = "bootstrapp"

process rasusaWithDocker {
    //container 'mbhall88/rasusa:latest'
    //containerOptions '--volume /Users/wasinee/cactus/Sporothrix_230320_NextSeq1000/Data/fastq:/data'

    output:
    path 'S01_sub.fastq'

    '''
        docker run -v ~/cactus/Sporothrix_230320_NextSeq1000/Data/fastq:/data  mbhall88/rasusa \
            /bin/rasusa --coverage 3 --genome-size 32Mb --input /data/S01_S1_R1_001.fastq.gz > S01_sub.fastq
    '''
}

process build_vlmc {
    publishDir "$params.bootstrapp/$sampleId"
    
    input:
    tuple val(rnd_seed), val(sampleId), path(bbduk_files)
    
    output:
    tuple val({rnd_seed}), path("*.bintree")

    """
        seqtk sample -s${rnd_seed} ${bbduk_files[0]} 4000000 > ${sampleId}_${rnd_seed}_sub1.fq
        seqtk sample -s${rnd_seed} ${bbduk_files[1]} 4000000 > ${sampleId}_${rnd_seed}_sub2.fq
        cat ${sampleId}_${rnd_seed}_sub1.fq ${sampleId}_${rnd_seed}_sub2.fq > ${sampleId}_merged.fastq
        dvstar \
        -m build \
        -p ${sampleId}_merged.fastq \
        --out-path ${sampleId}.bintree \
        --threshold $params.threshold --max-depth 4 --min-count 100 \
        --adjust-for-sequencing-errors \
        --sequencing-depth 20.0 \
        --sequencing-error-rate 0.001 \
        --pseudo-count-amount 0.01 \
        --temp-path tmp
        rm -rf ${sampleId}_${rnd_seed}_sub1.fq ${sampleId}_${rnd_seed}_sub2.fq ${sampleId}_merged.fastq
    """
}

rnd_ch = Channel.of(1000..1000000)
          .randomSample(1)
          .view()

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
        .view()

workflow {
    //rnd_ch.combine( reads_ch )
    //    .view()
    //| build_vlmc 
    //| view
    rasusaWithDocker | view
}