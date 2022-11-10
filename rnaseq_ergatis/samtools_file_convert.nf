#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

params.samtools_view_opts = ""
params.samtools_sort_opts = ""
//use --v for verbose summary
params.other_args = '--v'

process run_samtools_file_convert {
    outdir = "${workDir}/samtools_file_convert"
    publishDir outdir, mode: 'symlink', failOnError: true

    executor 'sge'
    memory '5 GB'

    input:
        path input_file
        val sortby
        val ref_fasta_file
        val samtools_bin_dir
    output:
        path "*.bam", glob: true, emit: bam_list
        path "*.sorted_by_position.bam", glob: true, emit: bam_sorted_by_position
        path "*.sorted_by_name.bam", glob: true, emit: bam_sorted_by_name
        path "*.sam", glob: true, emit: sam_list
        path "*.sorted_by_position.sam", glob: true, emit: sam_sorted_by_position
        path "*.sorted_by_name.sam", glob: true, emit: sam_sorted_by_name

    script:
        // string of options for file conversion (eg : 123)
        // 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM
        // "12" - sort by position, "13" - sort by name
        options = ""
        sort_opts = params.samtools_sort_opts
        if (sortby == "position") {
            options = "12"
        } else if (sortby == "name") {
            options = "13"
            if (!(params.samtools_sort_opts =~ /-n/)) {
                sort_opts = sort_opts +  ' -n'
            }
        }
        """
        /usr/bin/env perl ${params.bin_dir}/samtools_file_convert.pl \
            --infile=\$PWD/${input_file} \
            --infile_format=${input_file.extension}
            --options=${options}
            --reffile=\$PWD/${ref_fasta_file} \
            --outdir=\$PWD \
            --samtools_bin_dir=${samtools_bin_dir} \
            --samtools_view_options=${params.samtools_view_opts} \
            --samtools_sort_options=${sort_opts} \
            ${params.other_args}
        """
}

workflow samtools_file_convert {
    take:
        input_file
        sortby
        ref_fasta_file
        samtools_bin_dir
    main:
        run_samtools_file_convert(
            input_file
            , sortby
            , ref_fasta_file
            , samtools_bin_dir
            )
    emit:
        bam_list = run_samtools_file_convert.out.bam_list
        bam_sorted_by_position = run_samtools_file_convert.out.bam_sorted_by_position
        bam_sorted_by_name = run_samtools_file_convert.out.bam_sorted_by_name
        sam_list = run_samtools_file_convert.out.sam_list
        sam_sorted_by_position = run_samtools_file_convert.out.sam_sorted_by_position
        sam_sorted_by_name = run_samtools_file_convert.out.sam_sorted_by_name

}


// Workflow to call this as standalone
workflow {
    file_chunks_ch = channel.fromPath(params.fasta_file, checkIfExists:false, followLinks:true)
    if (! file_chunks_ch) {
        file_chunks_ch = channel.fromPath(params.fasta_list, checkIfExists:true, followLinks:true).splitText()
    }
    ref_fasta_file = channel.fromPath(params.ref_fasta_file, checkIfExists:false, followLinks:true)
    sortby = channel.fromValue(params.samtools_sortby)
    samtools_bin_dir = channel.fromPath(params.samtools_bin_dir, checkIfExists:true, followLinks:true)

    samtools_reference_index(
        file_chunks_ch
        , sortby
        , ref_fasta_file
        , samtools_bin_dir
        )
}