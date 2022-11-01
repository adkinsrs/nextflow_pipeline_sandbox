#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

//path to sample info file with information on all samples to be analyzed
params.sample_info = ""
//use --v for verbose summary
params.other_args = '--v'

process run_create_paired_list_file {
    outdir = "${workDir}/create_paired_list_file"
    publishDir outdir, mode: 'symlink', failOnError: true

    input:
        path list_file1
        path list_file2
    output:
        path "paired_input_file.list", glob: false

    """
    /usr/bin/env perl ${params.bin_dir}/create_paired_list_file.pl \
        --listfile1=\$PWD/${list_file1} \
        --listfile2=\$PWD/${list_file2} \
        --samplefile= ${params.sample_info} \
        --outdir=\$PWD \
        ${params.other_args}
    """
}

workflow create_paired_list_file {
    take: list_file1
    take: list_file2
    main:
        run_create_paired_list_file(
            list_file1
            , list_file2
            )
    emit: run_create_paired_list_file.out
}


// Workflow to call this as standalone
workflow {

    list_file1 = channel.fromPath(params.list_file1, checkIfExists:true, followLinks:true)
    list_file2 = channel.fromPath(params.list_file2, checkIfExists:true, followLinks:true)

    create_paired_list_file(
        list_file1
        , list_file2
        )
}