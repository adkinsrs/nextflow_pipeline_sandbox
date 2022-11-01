#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include {hisat2_build} from "./hisat2_build"
include {samtools_reference_index} from "./samtools_reference_index"
include {create_paired_list_file} from "./create_paired_list_file"

/*
process hisat2_build {
    input:
        path ref_fasta_file
        val hisat2_bin_dir
    output:
        path hisat2_build_index
        val prefix

}

process samtools_reference_index {
    input:
        path ref_fasta_file
        val samtools_bin_dir
    output:
        path samtools_index

}

process create_paired_list_file {
    input:
        path list_file1
        path list_file2
    output:
        path paired_list_file
}
*/

workflow {
    // In Ergatis, we could pass a file, a .list of files, or a directory.

    // this could just be read into their respective components directly, but it's easier to just do it here
    ref_fasta_file_ch = channel.fromPath(params.ref_fasta_file, checkIfExists:true, followLinks:true)

    list_file1_ch = channel.fromPath(params.list_file1, checkIfExists:true, followLinks:true)
    list_file2_ch = channel.fromPath(params.list_file2, checkIfExists:true, followLinks:true)

    hisat2_bin_dir = channel.fromPath(params.hisat2_bin_dir, checkIfExists:true, followLinks:true)
    samtools_bin_dir = channel.fromPath(params.samtools_bin_dir, checkIfExists:true, followLinks:true)

    main:
        hisat2_build(
            ref_fasta_file_ch
            , hisat2_bin_dir
            )
        samtools_reference_index(
            ref_fasta_file_ch
            , samtools_bin_dir
            )
        create_paired_list_file(
            list_file1_ch
            , list_file2_ch
        )
        fastqc_stats(create_paired_list_file.out)
        hisat2(
            create_paired_list_file.out
            , hisat2_build_index
            , hisat2_bin_dir
            , samtools_bin_dir
            )
        /*samtools_file_convert() // sorted by position
        samtools_file_convert() // sorted by name
        samtools_alignment_stats()
        align_hisat2_stats()
        percent_mapped_stats()
        bam2bigwig()
        rpkm_coverage_stats()
        wrapper_align()
        expression_plots()
        htseq()
        deseq()
        edgeR()
        filter_deseq()
        filter_edgeR()*/
    //emit: blah


}
