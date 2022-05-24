#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include workflow.config

params.abbreviation = "test"

process run_pmarks2bsml {

    input:
        path fasta_input
        val abbr

    output:
        path "${fasta_input.basename}.pmarks.bsml"

    """
    /usr/bin/env perl $params.bin_dir/pmarks2bsml \
        --fasta_input=${fasta_input} \
        --input_file=${fasta_input}.pmarks \
        --output=${fasta_input.basename}.pmarks.bsml \
        --project=${abbr} \
        --id_repository=${params.id_repository} \
    """
}

workflow pmarks2bsml {

    take: fasta_input, abbreviation

    run_pmarks2bsml(
        fasta_input
        , abbreviation
        )
    emit:
}

workflow {
    // In Ergatis, we could pass a file, a .list of files, or a directory.
    fasta_input_ch = channel.fromPath(params.fasta_file, checkIfExists=false, followLinks=true)
        .set{file_chunks_ch}
    fasta_list_ch = channel.fromPath(params.fasta_list, checkIfExists=false, followLinks=true)
        .splitText()
        .set{file_chunks_ch}

    abbreviation = channel.value(params.abbreviation)

    pmarks2bsml(file_chunks_ch, abbreviation)
}