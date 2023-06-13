#!/usr/bin/env nextflow

params.OUTPUT = "/Users/wasinee/nextflow/vlmc_distances/distance.tsv"

source = Channel.fromPath('./dvstar_output/S*.bintree')
                .ifEmpty { error "Cannot find any file matching" }
                .map { file -> tuple(file.baseName, file) }
                .toSortedList( { a, b -> a[0] <=> b[0] } ) // <=> is an operator for comparison
                
target = Channel.fromPath('./dvstar_output/S*.bintree')
                .ifEmpty { error "Cannot find any file matching" }
                .map { file -> tuple(file.baseName, file) }
                .toSortedList( { a, b -> a[0] <=> b[0] } ) // <=> is an operator for comparison
                

process foo {
  input:
    each source
    each target

  output:
    stdout

  shell:
  """
    #!/opt/anaconda3/envs/py310/bin/python
    VAL = 1
    id = "!{source[0]}"
    print("%d %s" % (VAL, id))
  """
}

process VLMCDistance {
    publishDir 'dvstar_distance'

    input:
        each source
        each target

    output:
    stdout

    shell:
    """
        #!/bin/bash
        
        DISTANCE=`dvstar \
                        --mode dissimilarity --dissimilarity dvstar \
                        --in-path !{source[1]} \
                        --to-path !{target[1]} \
                        | tail -n 1` 
        
       echo "!{source[0]}\t!{target[0]}\t\$DISTANCE" 2>/dev/null
    """
}

workflow {
    VLMCDistance(source, target) 
    | collectFile(name: 'dvstar_distance/distance.tsv', newLine: false)
}