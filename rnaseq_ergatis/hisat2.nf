#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line

//Scoring Options
// --mp MX,MN = Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. Default: MX = 6, MN = 2.
params.mismatch_penalties = "6,2"
//--sp MX,MN = Sets the maximum (MX) and minimum (MN) penalties for soft-clipping per base, both integers. Default: MX = 2, MN = 1.
params.softclip_penalties = "2,1"
//--rdg , = Sets the read gap open () and extend () penalties. A read gap of length N gets a penalty of + N * . Default: 5, 3.
params.read_gap_penalties = "5,3"
//--rfg , = Sets the reference gap open () and extend () penalties. A reference gap of length N gets a penalty of + N * . Default: 5, 3.
params.ref_gap_penalties = "5,3"
//--score-min = Sets a function governing the minimum alignment score needed for an alignment to be considered "valid". Default: L,0,-0.2.
params.min_alignment_score = "L,0,-0.2"
//Splice Alignment Options
// --pen-cansplice = Sets the penalty for each pair of canonical splice sites (e.g. GT/AG). Default: 0.
params.pen_can_splice = 0
//--pen-noncansplice = Sets the penalty for each pair of non-canonical splice sites (e.g. non-GT/AG). Default: 12.
params.pen_noncan_splice = 12
//--pen-canintronlen = Sets the penalty for long introns with canonical splice sites so that alignments with shorter introns are preferred to those with longer ones. Default: G,-8,1
params.pen_can_intron_len = "G,-8,1"
//--pen-noncanintronlen = Sets the penalty for long introns with noncanonical splice sites so that alignments with shorter introns are preferred to those with longer ones. Default: G,-8,1
params.pen_noncan_intron_len = "G,-8,1"
//The minimum intron length. Default:20.
params.min_intron_len = 20
//The maximum intron length. Default:500000.
params.max_intron_len = 500000
//--rna-strandness = Specify strand-specific information: the default is unstranded (leave field empty). F or R for single-end reads. FR or RF for paired-end reads.
params.rna_strandedness = "RF"
//Report alignments tailored specifically for Cufflinks.
params.dta_cufflinks = 1
//Reporting Options
// -k = number of distinct, primary alignments (default: 5 (HFM) or 10 (HGFM))
params.num_alignments = 5
//Paired-end Options
// -I/--minins = The minimum fragment length for valid paired-end alignments. Default: 0
params.min_fragment_len = 0
//-X/--maxins = The maximum fragment length for valid paired-end alignments. Default: 500
params.max_fragment_len = 500
//SAM options
// --no-unal = Suppress SAM records for reads that failed to align.
params.suppress_unalignments = 1
//Performance Options
// (-p) The number of threads to be used to align reads. Default:1.
params.num_threads = 1
// additional hisat2-build parameters
params.other_params = ''
//use --v for verbose summary
params.other_args = '--v'

process run_hisat2 {
    outdir = "${workDir}/hisat2"
    publishDir outdir, mode: 'symlink', failOnError: true

    executor 'sge'
    memory '5 GB'

    input:
        tuple val(seq_file1), val(seq_file2)
        path hisat2_build_index
        val hisat2_build_prefix
        val hisat2_bin_dir
        val samtools_bin_dir

    output:
        path "*.accepted_hits.bam",  glob: true, emit: bam_files
        path "*.accepted_hits.sam",  glob: true, emit: sam_files

    """
    /usr/bin/env perl ${params.bin_dir}/hisat2.pl \
        --seq1file=${seq_file1} \
        --seq2file=${seq_file2} \
        --hisat2_index_dir=${hisat2_build_index} \
        --prefix=${hisat2_build_prefix} \
        --outdir=\$PWD \
        --mismatch-penalties=${params.mismatch_penalties} \
        --softclip-penalties=${params.softclip_penalties} \
        --read-gap-penalties=${params.read_gap_penalties} \
        --ref-gap-penalties=${params.ref_gap_penalties} \
        --min-intronlen=${params.min_intron_len} \
        --max-intronlen=${params.max_intron_len} \
        --score-min=${params.min_alignment_score} \
        --pen-cansplice=${params.pen_can_splice} \
        --pen-noncansplice=${params.pen_noncan_splice} \
        --pen-canintronlen=${params.pen_can_intron_len} \
        --pen-noncanintronlen=${params.pen_noncan_intron_len} \
        --rna-strandness=${params.rna_strandedness} \
        --dta-cufflinks=${params.dta_cufflinks} \
        --num-alignments=${params.num_alignments} \
        --minins=${params.min_fragment_len} \
        --maxins=${params.max_fragment_len} \
        --no-unal=${params.suppress_unalignments} \
        --num-threads=${params.num_threads} \
        --hisat2_bin_dir=${hisat2_bin_dir} \
        --samtools_bin_dir=${samtools_bin_dir} \
        --args=${params.other_params} \
        ${params.other_args}
    """
}

workflow hisat2 {
    take:
        paired_files
        hisat2_build_index
        hisat2_build_prefix
        hisat2_bin_dir
        samtools_bin_dir
    main:
        paired_files
            .splitText()
            .splitCsv(sep:"\t", header:false)
            .map { row -> tuple(row[0], row[1]) }
        run_hisat2(
            paired_files
                .splitText()
                .splitCsv(sep:"\t", header:false)
                .map { row -> tuple(row[0], row[1]) }
            , hisat2_build_index
            , hisat2_build_prefix
            , hisat2_bin_dir
            , samtools_bin_dir
            )
    emit:
       bam_files = run_hisat2.out.bam_files
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

    hisat2_index_ch = channel.fromPath(params.hisat2_index_dir, checkIfExists:true, followLinks:true)
    hisat2_index_prefix_ch = channel.value(params.hisat2_prefix, checkIfExists:true, followLinks:true)

    hisat2_bin_dir = channel.fromPath(params.hisat2_bin_dir, checkIfExists:true, followLinks:true)
    samtools_bin_dir = channel.fromPath(params.samtools_bin_dir, checkIfExists:true, followLinks:true)

    hisat2 (
        paired_input_file_ch
        , hisat2_index_ch
        , hisat2_index_prefix_ch
        , hisat2_bin_dir
        , samtools_bin_dir
        )
}
