"""
This part of the workflow handles sorting downloaded sequences and metadata
by aligning them to reference sequences.

It produces output files as

    metadata = "data/metadata.tsv"
    sequences = "data/sequences.fasta"

"""

rule sort:
    input:
        sequences = rules.curate.output.sequences
    output:
        "results/sequences.fasta",
    shell:
        '''
        seqkit rmdup {input.sequences} > {output}
        '''

rule metadata:
    input:
        # via add_metadata_columns so the `url` column is present
        metadata = rules.add_metadata_columns.output.metadata,
        sequences = "results/sequences.fasta"
    output:
        metadata = "results/metadata_raw.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade:
    input:
        sequences = "results/sequences.fasta",
        ref = "data/references/reference.fasta"
    output:
        nextclade = "results/nextclade.tsv"
    params:
        dataset = "data/references/",
        output_columns = "seqName clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic",
        min_length = config["curate"]["min_length"]
    log:
        "logs/nextclade.txt"
    threads: workflow.cores
    shell:
        """
        nextclade3 run -D {params.dataset}  -j {threads} \
                          --output-columns-selection {params.output_columns} \
                          --min-length {params.min_length} \
                          --output-tsv {output.nextclade} \
                          {input.sequences} > {log} 2>&1
        """

rule extend_metadata: 
    input:
        nextclade = "results/nextclade.tsv",
        metadata = "results/metadata_raw.tsv"
    output:
        metadata = "results/metadata.tsv"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """

