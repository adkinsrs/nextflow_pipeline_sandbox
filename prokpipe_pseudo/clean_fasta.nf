#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// These can be overwritten by passing --<param_name> in the command line
params.clean_fasta_opts = "-u"
params.nucleotide = true
nuc = params.nucleotide ? "T" : "F"

process remove_cr_line_breaks {
    outdir = "$workDir/clean_fasta"

    publishDir outdir, mode: 'symlink', failOnError: true

    input:
        path input_file // should be absolute values
    output:
        // simpleName - https://www.nextflow.io/docs/latest/script.html#check-file-attributes
        // Note will not work if file basename contains periods
        path "${input_file.simpleName}.tmp.fa", glob: false

    """
    /usr/bin/env perl -p -w -e "s/\r//g" $input_file > "${input_file.simpleName}.tmp.fa"
    """
}

process run_clean_fasta {
    outdir = "$workDir/clean_fasta"

    publishDir outdir, mode: 'symlink', failOnError: true

    input:
        path fasta_input

    output:
        path "${fasta_input.simpleName}.fa", glob: false

    """
    /usr/bin/env perl ${params.bin_dir}/clean_fasta.pl \
        ${params.clean_fasta_opts} \
        $fasta_input \
        --output=${fasta_input.simpleName}.fa
    """
}

process replace_gaps_in_nucleotides {
    outdir = "${workDir}/final_clean"

    publishDir outdir, mode: 'symlink', failOnError: true

    input:
        path fasta_input

    output:
        path "${fasta_input.simpleName}.cleaned.fa", glob: false

    """
    /usr/bin/env perl ${params.bin_dir}/replace_nuc_gaps.pl \
        --nucleotide=$nuc \
        --input_file=$fasta_input \
        --output_file=${fasta_input.simpleName}.cleaned.fa
    """
}

workflow clean_fasta {
    take: fasta_input
    main:
        remove_cr_line_breaks(fasta_input)
        run_clean_fasta(remove_cr_line_breaks.out)
        replace_gaps_in_nucleotides(run_clean_fasta.out)
    emit: replace_gaps_in_nucleotides.out
}


// Workflow to call this as standalone
workflow {
    // In Ergatis, we could pass a file, a .list of files, or a directory.
    file_chunks_ch = channel.fromPath(params.fasta_file, checkIfExists:false, followLinks:true)
    if (! file_chunks_ch) {
        file_chunks_ch = channel.fromPath(params.fasta_list, checkIfExists:true, followLinks:true).splitText()
    }

    clean_fasta(
        file_chunks_ch
        )
}