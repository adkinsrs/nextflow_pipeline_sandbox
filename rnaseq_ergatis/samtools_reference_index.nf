#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

//use --v for verbose summary
params.other_args = '--v'

process run_samtools_reference_index {
    outdir = "${workDir}/samtools_reference_index"
    publishDir outdir, mode: 'symlink', failOnError: true

    input:
        path ref_fasta_file
        path samtools_bin_dir
    output:
        path "${ref_fasta_file.simpleName}.fa.fai", glob: false


    """
    /usr/bin/env perl ${params.bin_dir}/samtools_reference_index.pl \
        --reffile=${ref_fasta_file} \
        --outdir=${outdir} \
        --samtools_bin_dir=${samtools_bin_dir} \
        ${params.other_args}
    """
}

workflow samtools_reference_index {
    take: ref_fasta_file, samtools_bin_dir
    main:
        run_samtools_reference_index(
            ref_fasta_file
            , samtools_bin_dir
            )
    emit: run_samtools_reference_index.out
}


// Workflow to call this as standalone
workflow {
    ref_fasta_file = channel.fromPath(params.ref_fasta_file, checkIfExists:true, followLinks:true)
    samtools_bin_dir = channel.fromPath(params.samtools_bin_dir, checkIfExists:true, followLinks:true)


    samtools_reference_index(
        ref_fasta_file
        , samtools_bin_dir
        )
}