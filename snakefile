###############
# Snakemake execution templates:

# To run a default protein xy run:
# snakemake  auspice/<your_virus>_protein_xy.json --cores 9

# To run a default whole genome run (>6400bp):
# snakemake auspice/<your_virus>_genome.json --cores 9

###############
wildcard_constraints:
    seg="protein_xy|genome"  # Define segments to analyze, e.g. vp1, whole-genome. This wildcard will be used in the rules "{seg}" to define the path or protein to use
   
# Define segments to analyze
segments = ['protein_xy', 'genome'] # This is only for the expand in rule all

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/<your_virus>_{segs}.json", segs=segments) ## TODO: replace <your_virus> with actual virus name (Ctrl+H)


# Rule to handle configuration files and data file paths
rule files:
    input:
        sequence_length =   "{seg}",
        colors =            "config/colors.tsv",
        dropped_strains =   "config/dropped_strains.txt",
        regions=            "config/geo_regions.tsv",
        lat_longs =         "config/lat_longs.tsv",
        reference =         "{seg}/config/reference_sequence.gb", ####TODO: provide a reference sequence
        auspice_config =    "{seg}/config/auspice_config.json",
        clades =            "{seg}/config/clades_genome.tsv",
        meta=               "data/metadata.tsv",
        extended_metafile=  "data/meta_manual.tsv", ###TODO: Add an empty tsv file to this path or metadata for one of your sequences


files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

rule fetch:  ####TODO: go to ingest readme and read through instructions and files needed!
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata=files.meta
    params:
        seq="ingest/data/sequences.fasta",
        meta="ingest/data/metadata.tsv"
    shell:
        """
        cd {input.dir} 
        snakemake --cores 9 all
        cd ../
        cp -u {params.seq} {output.sequences}
        cp -u {params.meta} {output.metadata}
        """

##############################
# Update strain names
# If strain name == accession -> fetching real strain names from genbank
# Depending on how many sequences you have, it will run for a long time! >30min. Comment out to skip!
###############################

rule update_strain_names:
    message:
        """
        Updating strain information in metadata.
        """
    input:
        file_in =  files.meta
    params:
        backup = "data/strain_names_previous_run.tsv" ####TODO: provide an empty file for first run
    output:
        file_out = "data/updated_strain_names.tsv"
    shell:
        """
        time bash scripts/update_strain.sh {input.file_in} {params.backup} {output.file_out}
        """


##############################
# BLAST
# blast fasta files for your specific proteins
# cut out your protein from fasta sequences
###############################

rule blast:
    input: 
        blast_db_file = "data/references/reference_blast.fasta",    ####TODO: provide a blast reference
        seqs_to_blast = rules.update_sequences.output.sequences
    output:
        blast_out = "temp/blast_out.csv"
    params:
        blast_db = "temp/entero_db_vp1" ##TODO: rename
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort: #TODO: change the parameters in blast_sort.py (replace vp1 with your specific protein)
    input:
        blast_result = rules.blast.output.blast_out, # output blast (for your protein)
        input_seqs = rules.update_sequences.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        protein = [600,915], #TODO: min & max length for protein
        whole_genome = [6400,8000], #TODO: min & max length for whole genome
        range = "{seg}" # this is determining the path it takes in blast_sort (protein-specific or whole genome)
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --protein_length {params.protein}  --whole_genome_length {params.whole_genome} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range}
        """
##############################
# AUGUR CURATE AND MERGE
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
# Merge with other metadata files you might have
###############################

rule curate_meta_dates:
    message:
        """
        Cleaning up metadata with augur curate and merge your metadata with the one from ingest
        """
    input:
        meta = files.meta,
        metadata=files.extended_metafile,  ###TODO: Add an empty tsv file to this path or metadata for one of your sequence.
    params:
        strain_id_field="accession",
        date_column="date",
        format=['%Y', '%m.%Y', '%d.%m.%Y', "%b-%Y", "%d-%b-%Y","%Y-%m-%d"],
        temp_metadata="data/temp_curated.tsv"  # Temporary file
    output:
        meta = "data/metadata_ingest_curated.tsv",
        metadata = "data/metadata_manual_curated.tsv",
        final_metadata="data/metadata_curated_merged.tsv"  # Final merged output file

    shell:
        """
        # Normalize strings for metadata
        augur curate normalize-strings --metadata {input.metadata} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata}

        # Format dates for metadata
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --no-mask-failure \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        
        # Remove temporary file
        rm {params.temp_metadata}

        # Normalize strings for genbank metadata
        augur curate normalize-strings --metadata {input.meta} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata}

        # Format dates for genbank metadata
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --no-mask-failure \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.meta}
        
        # Remove temporary file
        rm {params.temp_metadata}

        # Merge curated metadata
        augur merge --metadata meta={output.meta} extended_meta={output.metadata}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.final_metadata}
        """

##############################
# Indexing sequences and filter them.
###############################


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.blast_sort.output.sequences
    output:
        sequence_index = "{seg}/results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata =rules.curate_meta_dates.output.final_metadata,
        exclude = files.dropped_strains
    output:
        sequences = "{seg}/results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 4000, # add a limit per group
        strain_id_field= "accession",
        min_date = 1950  # add a reasonable min date
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --output {output.sequences}
        """
# --exclude-where ... or other parameters can be added, see `augur filter --h` for more options

##############################
# Reference for alignment added to sub-folders
###############################

rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format and putting it in the reference folders of your proteins
        """
    input:
        reference = files.reference

    output:
        reference = "{seg}/results/reference_sequence.fasta"
    run:
        from Bio import SeqIO 
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign.
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "{seg}/results/aligned.fasta"

    params:
            nuc_mismatch_all = 10,
            nuc_seed_length = 30
    shell:
        """
        nextclade run \
        {input.sequences}  \
        --input-ref {input.reference}\
        --allowed-mismatches {params.nuc_mismatch_all} \
        --min-length {params.nuc_seed_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
        """

##############################
# Building a tree
###############################

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.align.output.alignment

    output:
        tree = "{seg}/results/tree_raw.nwk"

    threads: 9
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree}
        """

##############################
# Refine to a timeline
###############################

rule refine:
    message:
        """
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.curate_meta_dates.output.final_metadata,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 3, # set to 6 if you want more control over outliers
        strain_id_field ="accession",
        clock_rate = 0.004, # remove for estimation by augur; check literature
        clock_std_dev = 0.0015

    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --clock-rate {params.clock_rate}\
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

##############################
# Ancestral sequences and amino acids
###############################

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment

    output:
        node_data = "{seg}/results/nt_muts.json"

    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous\
            --inference {params.inference}
        """
 
rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "{seg}/results/aa_muts.json"

    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata =rules.curate_meta_dates.output.final_metadata
    output:
        node_data = "{seg}/results/traits.json"
        
    params:
        traits = "country",
        strain_id_field= "accession"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.traits} \
            --confidence
        """

##############################
# Assign clades or subgenotypes based on list provided
###############################
rule clades: 
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades # TODO: assign mutations to specific clades
    output:
        clade_data = "{seg}/results/clades.json"

    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

#########################
#  EXPORT
#########################
rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.curate_meta_dates.output.final_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= "accession"

    output:
        auspice_json = "auspice/<your_virus>_{seg}-accession.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
                {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
        
###############################
# Change from accession to strain name view in tree
################################

rule rename_json:
    input:
        auspice_json= rules.export.output.auspice_json,
        metadata =rules.curate_meta_dates.output.final_metadata,
    output:
        auspice_json="auspice/<your_virus>_{seg}.json"
    params:
        strain_id_field="accession",
        display_strain_field= "strain"
    shell:
        """
        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --metadata-id-columns {params.strain_id_field} \
                --input-auspice-json {input.auspice_json} \
                --display-strain-name {params.display_strain_field} \
                --output {output.auspice_json}

        mkdir -p auspice/accession/ && mv {input.auspice_json} auspice/accession/
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results/*",
        # "auspice/*",
        "temp/*"
    shell:
        "rm -rfv {params}"