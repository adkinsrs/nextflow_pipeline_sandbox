#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

// additional hisat2-build parameters
params.other_params = ''
//use --v for verbose summary
params.other_args = '--v'

process run_hisat2_build {
    outdir = "${workDir}/hisat2_build"
    publishDir outdir, mode: 'symlink', failOnError: true

    executor 'sge'
    memory '5 GB'

    input:
        path ref_fasta_file
        val hisat2_bin_dir
    output:
        path "index", emit: index   // Need to pass the directory into subsequent components
        val ref_fasta_file.baseName, emit: prefix

    // script makes output directory if it does not exist
    // Script expects full path for ref fasta file
    """
    mkdir \$PWD/index
    /usr/bin/env perl ${params.bin_dir}/hisat2_build.pl \
        --reffile=\$PWD/${ref_fasta_file} \
        --outdir=\$PWD/index \
        --prefix=${ref_fasta_file.baseName} \
        --hisat2_bin_dir=${hisat2_bin_dir} \
        --args=${params.other_params} \
        ${params.other_args}
    """
}

workflow hisat2_build {
    // Only one input channel per "take" and one output channel per "emit"
    take:
        ref_fasta_file
        hisat2_bin_dir
    main:
        run_hisat2_build(
            ref_fasta_file
            , hisat2_bin_dir
            )
    emit:
        index = run_hisat2_build.out.index
        prefix = run_hisat2_build.out.prefix
}


// Workflow to call this as standalone
workflow {
    ref_fasta_file = channel.fromPath(params.ref_fasta_file, checkIfExists:true, followLinks:true)
    hisat2_bin_dir = channel.fromPath(params.hisat2_bin_dir, checkIfExists:true, followLinks:true)

    hisat2_build(
        ref_fasta_file
        , hisat2_bin_dir
        )
}