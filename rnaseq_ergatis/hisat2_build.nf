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

    input:
        path ref_fasta_file
        path hisat2_bin_dir
    output:
        path "index", emit: index   // Need to pass the directory into subsequent components
        val ref_fasta_file.simpleName, emit: prefix



    """
    mkdir ${outdir}/hisat2
    /usr/bin/env perl ${params.bin_dir}/hisat2_build.pl \
        --reffile=${ref_fasta_file} \
        --outdir=${outdir}/index \
        --prefix=${ref_fasta_file.simpleName} \
        --hisat2_bin_dir=${hisat2_bin_dir} \
        --args=${params.other_params} \
        ${params.other_args}
    """
}

workflow hisat2_build {
    take: ref_fasta_file, hisat2_bin_dir
    main:
        run_hisat2_build(
            ref_fasta_file
            , hisat2_bin_dir
            )
    emit: run_hisat2_build.out.index, run_hisat2_build.out.prefix
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