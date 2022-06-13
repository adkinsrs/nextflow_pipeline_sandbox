#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

includeConfig workflow.config

include {clean_fasta} from "./clean_fasta"
include {create_pseudomolecules} from "./create_pseudomolecules"
include {pmarks2bsml} from "./pmarks2bsml"

params.abbreviation = "test"    // Codename to use in output files

process clean_fasta {
    input:
        path file_chunks_ch
    output:
        path clean_fasta_out

}

process create_pseudomolecules {

}

process pmarks2bsml {

}

workflow {
    // In Ergatis, we could pass a file, a .list of files, or a directory.
    file_chunks_ch = channel.fromPath(params.fasta_file, checkIfExists:false, followLinks:true)
    if (! file_chunks_ch) {
        file_chunks_ch = channel.fromPath(params.fasta_list, checkIfExists:true, followLinks:true).splitText()
    }

    main:
        // Executes per file
        clean_fasta(
            file_chunks_ch
            )
        create_pseudomolecules(
            clean_fasta.out
            , params.accession_file
            , params.abbreviation
            )
        pmarks2bsml(
            create_pseudomolecules.out
            , params.abbreviation
            )
    emit: pmarks2bsml.out

}
