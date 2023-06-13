#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//adapters_file = file(params.adapters)
dataDir = "/Users/wasinee/cactus/Sporothrix_230320_NextSeq1000/Data/fastq" 
params.reads = "/Users/wasinee/cactus/Sporothrix_230320_NextSeq1000/Data/fastq/S02_*{R1,R2}*.gz"
params.postTrim = "postTrim"

println "reads: $params.reads"

//Creates the `read_pairs` channel that emits for each read-pair a tuple containing
//three elements: the pair ID, the first read-pair file and the second read-pair file

process foo {
  debug true

  input:
  tuple val(sampleId), file(reads)

  script:
  """
  echo your_command --sample $sampleId --reads $reads
  """
}

//process 1: Adapter and quality trimming

process bbduk {
    tag {sampleId}
    publishDir "$params.postTrim/$sampleId"
    debug true

    input:
    tuple val(sampleId), file(reads)

    output:
    tuple val({sampleId}), path("*.bbduk.fastq")
    
    script:
    """
        bbduk.sh \
             in1=${reads[0]} \
             in2=${reads[1]} \
             ref="adapters" \
             out1=${reads[0].baseName}.bbduk.fastq \
             out2=${reads[1].baseName}.bbduk.fastq \
             stats=${sampleId}.stats.txt \
             threads=${task.cpus} \
             minlen=30 \
             qtrim=rl \
             trimq=10 \
             ktrim=r \
             k=30 \
             mink=11 \
             hdist=1 \
             trimbyoverlap \
             trimpairsevenly \
    """
}

process FastQC {
    input:
    path bbduk_files

    script:
    """
        mkdir fastqc_output
        fastqc -o fastqc_output bbduk_files
    """
}

process Spades {
        publishDir 'spades_output'

        input:
        path bbduk_files

        output:
        path "spades_output/*.{fasta,fastq}"


        """
        spades.py \
        -k 21,33,55,77,99,127 \
        --careful \
        --threads 4 \
        --pe1-1 ${bbduk_files[0]} \
        --pe1-2 ${bbduk_files[1]} \
        -o spades_output \
        """
}

process merge_fasta {
    publishDir "$params.postTrim"

    input:
    tuple val(sampleId), path(bbduk_files)

    output:
    path "*.fastq"

    """
        cat ${bbduk_files[0]} ${bbduk_files[1]} > ${sampleId}_merged.fastq
    """
}

workflow {
    post_trim = Channel
        .fromFilePairs('/Users/wasinee/cactus/Sporothrix_230320_NextSeq1000/Data/fastq/S*{R1,R2}*.gz')
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
        | bbduk | view
}