#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

// path to fastqc package binary directory
params.fastqc_bin_dir = "/usr/local/packages/fastqc/"

params.other_params = ''

process run_fastqc_stats {
    outdir = "${workDir}/fastqc_stats"
    publishDir outdir, mode: 'symlink', failOnError: true

    input:
        tuple file(seq_file1), file(seq_file2)
    output:
        path "*.per_base_quality.png", glob: true, emit: per_base_quality
        path "*.sequence_length_distribution.png", glob: true, emit: sequence_length_distribution

    """
    /usr/bin/env perl ${params.bin_dir}/fastqc_stats.pl \
        --seq1file=${seq_file1} \
        --seq2file=${seq_file2} \
        --outdir=${outdir} \
        --fastqc_bin_dir=${params.fastq_bin_dir}
        --args=${params.other_params}
    """
}

workflow fastqc_stats {
    take: paired_files
    main:
        run_fastqc_stats (
            paired_files
            )
    emit: run_fastqc_stats.out.per_base_quality, run_fastqc_stats.out.sequence_length_distribution
}


// Workflow to call this as standalone
workflow {

    // In this component for Ergatis, we could pass a pair of individual files or a "paired files" list
    list_file1_ch = channel.fromPath(params.list_file1, checkIfExists:true, followLinks:true)
    list_file2_ch = channel.fromPath(params.list_file2, checkIfExists:true, followLinks:true)
    paired_input_file_ch = channel.fromFilePairs()
    if (! (list_file1_ch && list_file2_ch)) {
        paired_input_file_ch = channel.fromPath(params.paired_input_file_list, checkIfExists:true, followLinks:true).splitText()
    }


    fastqc_stats (
            paired_input_file_ch
        )
}