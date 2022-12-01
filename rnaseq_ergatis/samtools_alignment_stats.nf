#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

//use --v for verbose summary
params.other_args = '--v'

process run_samtools_alignment_stats {
    outdir = "${workDir}/samtools_reference_index"
    publishDir outdir, mode: 'symlink', failOnError: true

    executor 'sge'
    memory '5 GB'

    input:
        path bam_file
        val samtools_bin_dir
    output:
        path "*.mapstats.txt", glob: true

    """
    /usr/bin/env perl ${params.bin_dir}/samtools_alignment_stats.pl \
        --infile=\$PWD/${bam_file} \
        --outdir=\$PWD \
        --samtools_bin_dir=${samtools_bin_dir} \
        ${params.other_args}
    """
}

workflow samtools_alignment_stats {
    take:
        bam_file
        samtools_bin_dir
    main:
        run_samtools_alignment_stats(
            bam_file
            , samtools_bin_dir
            )
    emit: run_samtools_alignment_stats.out
}


// Workflow to call this as standalone
workflow {
    file_chunks_ch = channel.fromPath(params.fasta_file, checkIfExists:false, followLinks:true)
    if (! file_chunks_ch) {
        file_chunks_ch = channel.fromPath(params.fasta_list, checkIfExists:true, followLinks:true).splitText()
    }
    samtools_bin_dir = channel.fromPath(params.samtools_bin_dir, checkIfExists:true, followLinks:true)

    samtools_reference_index(
        file_chunks_ch
        , samtools_bin_dir
        )
}