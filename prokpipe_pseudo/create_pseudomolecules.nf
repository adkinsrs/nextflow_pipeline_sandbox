#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include workflow.config


params.contig_input = "/path/to/contigs.fasta"
params.accession_input = "/path/to/accessions.txt"
params.abbreviation = "test"

params.format = "fasta"
params.nucmer_config = "-c 100 -maxmatch "
params.database = "nucleotide"
params.linker_sequences = null
params.debug=4

process run_create_pseudomolecules {
    scratch true

    input:
        path contig_input
        path accession_input
        val abbr

    output:
        path output_dir


    """
    /usr/bin/env perl $params.bin_dir/create_pseudomolecules.pl \
        --input_file=$accession_input \
        --contig_input=$contig_input \
        --strain=$abbr \
        --database=$params.database \
        --format=$params.format \
        --output_dir=$output_dir \
        --nucmer_exec=$params.nucmer_exec \
        --coords_exec=$params.coords_exec \
        --config_param=$params.nucmer_config \
        --linker_sequences=$params.linker_sequences \
        --log=$output_dir/create_pseudomolecules.log \
        --debug=$params.debug \
    """
}


workflow create_pseudomolecules {
    take: contig_input, accession_input, abbreviation
    run_create_pseudomolecules(
        contig_input
        , accession_input
        , abbreviation
        )
    emit: run_create_pseudomolecules.out

}

workflow {
    contig_input = channel.fromPath(params.contig_input)
    accession_input = channel.fromPath(params.accession_input)
    // SAdkins - not sure if channel.value is necessary vs just using params.abbreviation only
    abbreviation = channel.value(params.abbreviation)

    create_pseudomolecules(
        contig_input
        , accession_input
        , abbreviation
        )
}