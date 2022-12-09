#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line
params.other_params = ''

process run_fastqc_stats {
    outdir = "${workDir}/fastqc_stats"
    publishDir outdir, mode: 'symlink', failOnError: true

    executor 'sge'
    memory '5 GB'

    input:
        // comma-separated paths are perfectly valid input, which means using "path" types to stage is not possible
        tuple val(seq_file1), val(seq_file2)
        each fastqc_bin_dir
    output:
        path "*.per_base_quality.png", glob: true, emit: per_base_quality
        path "*.sequence_length_distribution.png", glob: true, emit: sequence_length_distribution

    """
    /usr/bin/env perl ${params.bin_dir}/fastqc_stats.pl \
        --seq1file=\$PWD/${seq_file1} \
        --seq2file=\$PWD/${seq_file2} \
        --outdir=\$PWD \
        --fastqc_bin_dir=${fastqc_bin_dir}
        --args=${params.other_params}
    """
}

workflow fastqc_stats {
    take:
        paired_files
        fastqc_bin_dir
    main:
        run_fastqc_stats (
            paired_files
                .splitText()
                .splitCsv(sep:"\t", header:false)
                .map { row -> tuple(row[0], row[1]) }   // can instead use { tuple(it[0], it[1]) } as "it" is implicit
            , fastqc_bin_dir
            )
    emit:
        base_quality = run_fastqc_stats.out.per_base_quality.collect()
        seq_len_distribution = run_fastqc_stats.out.sequence_length_distribution.collect()
}


// Workflow to call this as standalone
workflow {
    // In this component for Ergatis, we could pass a pair of individual files or a "paired files" list
    paired_input_file_ch = channel.fromPath(params.paired_input_file_list, checkIfExists:true, followLinks:true)
    if (! paired_input_file_ch) {
        //TODO: zip, tab-delimit, and process each pair
        list_file1_ch = channel.fromPath(params.list_file1, checkIfExists:true, followLinks:true).splitText()
        list_file2_ch = channel.fromPath(params.list_file2, checkIfExists:true, followLinks:true).splitText()
        paired_input_file_ch = list_file1_ch.merge(list_file2_ch)
    }

    fastqc_stats (
            paired_input_file_ch
            , fastqc_bin_dir
        )
}
