#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include workflow.config

// These can be overwritten by passing --<param_name> in the command line
params = {
    clean_fasta_opts = "-u"
    nucleotide = true
    outdir = file("$workDir/clean_fasta")
}

process remove_cr_line_breaks {

    params.outdir.mkdir()

    input:
        path input_file // should be absolute values
    output:
        path "$outdir/${input_file.basename}.tmp.fa"

    """
    /usr/bin/env perl -p -w -e "s/\r//g" $input_file; > "${params.outdir}/${input_file.basename}.tmp.fa"
    """

}

process run_clean_fasta {
    scratch true

    input:
        path fasta_input

    output:
        path "${params.outdir}${fasta_input.basename}.clean.fa"

    """
    /usr/bin/env perl $params.bin_dir/clean_fasta.pl \
    $params.clean_fasta_opts \
    $fasta_input \
    --output=${params.outdir}/${fasta_input.basename}.clean.fa" \

    """

}

process replace_gaps_in_nucleotides {
    outdir = file("$workDir/final_clean").mkdir()

    scratch true

    nuc="F"
    if (params.nucleotide) {
        nuc="T"
    }

    input:
        path fasta_input

    output:
        path "${outdir}/${fasta_input.basename}.fa"

    """
    /usr/bin/env perl $params.bin_dir/replace_nuc_gaps.pl \
        --nucleotide=$nuc \
        --input_file=$fasta_input \
        --output_file=${outdir}/${fasta_input.basename}.fa" \
    """

}

workflow clean_fasta {
    take: fasta_input
    main:
        remove_cr_line_breaks(fasta_input)
        run_clean_fasta(remove_cr_line_breaks.out)
        replace_gaps_in_nucleotides(clean_fasta.out)
    emit: replace_gaps_in_nucleotides.out
}


// Workflow to call this as standalone
workflow {
    // In Ergatis, we could pass a file, a .list of files, or a directory.
    fasta_input_ch = channel.fromPath(params.fasta_file, checkIfExists=false, followLinks=true)
        .set{file_chunks_ch}
    fasta_list_ch = channel.fromPath(params.fasta_list, checkIfExists=false, followLinks=true)
        .splitText()
        .set{file_chunks_ch}

    clean_fasta(
        file_chunks_ch
        )
}