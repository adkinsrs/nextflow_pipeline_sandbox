#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include {hisat2_build} from "./hisat2_build"
include {samtools_reference_index} from "./samtools_reference_index"
include {create_paired_list_file} from "./create_paired_list_file"
include {fastqc_stats} from "./fastqc_stats"
include {hisat2} from "./hisat2"
include {samtools_file_convert} from "./samtools_file_convert"

workflow samtools_file_convert_position {
    take:
        input_file
        samtools_bin_dir
    main:
        samtools_file_convert(
            input_file
            , "position"
            , ""
            , samtools_bin_dir
            )   // sorted by position
    emit:
        bam_sorted_by_position = samtools_file_convert.out.bam_sorted_by_position

}
workflow samtools_file_convert_name {
    take:
        input_file
        samtools_bin_dir
    main:
        samtools_file_convert(
            input_file
            , "name"
            , ""
            , samtools_bin_dir
            )   // sorted by name
    emit:
        bam_sorted_by_name = samtools_file_convert.out.bam_sorted_by_name
}

workflow {
    // In Ergatis, we could pass a file, a .list of files, or a directory.

    // this could just be read into their respective components directly, but it's easier to just do it here
    ref_fasta_file_ch = channel.fromPath(params.ref_fasta_file, checkIfExists:true, followLinks:true)

    list_file1_ch = channel.fromPath(params.list_file1, checkIfExists:true, followLinks:true)
    list_file2_ch = channel.fromPath(params.list_file2, checkIfExists:true, followLinks:true)

    hisat2_bin_dir = channel.fromPath(params.hisat2_bin_dir, checkIfExists:true, followLinks:true)
    samtools_bin_dir = channel.fromPath(params.samtools_bin_dir, checkIfExists:true, followLinks:true)
    fastqc_bin_dir = channel.fromPath(params.fastqc_bin_dir, checkIfExists:true, followLinks:true)

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
        /* NOTE: This tool fails sometimes
        fastqc_stats(
            create_paired_list_file.out
            , fastqc_bin_dir
            )
        */
        hisat2(
            create_paired_list_file.out
            , hisat2_build.out.index
            , hisat2_build.out.prefix
            , hisat2_bin_dir
            , samtools_bin_dir
            )

        samtools_file_convert_position(
            hisat2.out.bam_files
            , samtools_bin_dir
            )   // sorted by position
        samtools_file_convert_name(
            hisat2.out.bam_files
            , samtools_bin_dir
            )   // sorted by name
        samtools_alignment_stats(
            samtools_file_convert_position.out
            , samtools_bin_dir
        )
        /*
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
