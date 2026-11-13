###############
# Snakemake execution templates:

# To run a default protein xy run:
# snakemake  auspice/YOUR_VIRUS_protein_xy.json --cores 9

# To run a default whole genome run (>6400bp):
# snakemake auspice/YOUR_VIRUS_genome.json --cores 9

import os
from datetime import date

TAXID = "<TAXID>"  # NCBI Taxonomy ID for your virus


# Load config file
if not config:
    configfile: "config/config.yaml"

# Load environment variables
# Try to load .env, but don't fail if it doesn't exist (for Actions)
try:
    from dotenv import load_dotenv
    load_dotenv(".env")
except:
    pass


###############
#ensure protein_xy name similar to that found in the reference_sequence.gb CDS
wildcard_constraints:
    seg="protein_xy|protein_yz|genome"  # Define segments to analyze, e.g. vp1, whole-genome. This wildcard will be used in the rules "{seg}" to define the path or protein to use

# Define segments to analyze
segments = ['protein_xy', 'protein_yz','genome'] # This is only for the expand in rule all. TODO: replace <protein_xy> with "vp1" or another protein for which you would like to have a separate workflow

# parameters
FETCH_SEQUENCES = True
CODING_GENES = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"] # list all coding genes - it is used for augur ancestral to reconstruct the ancestral sequences for each gene. 

# Rule to handle configuration files
rule files:
    input:
        sequence_length =   "{seg}",
        exclude =           "config/exclude.txt",
        include =           "config/include.txt",
        reference =         "config/reference.gbk",
        gff_reference =     "{seg}/config/annotation.gff3",
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv",
        meta_collab =       "data/meta_collab.tsv",
        SEQUENCES =         "data/fetch/sequences.fasta",
        METADATA =          "data/fetch/metadata.tsv",

files = rules.files.input
##############################

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/YOUR_VIRUS_{segs}.json", segs=segments), ## TODO: replace YOUR_VIRUS with actual virus name (Ctrl+H). Typical naming convention: e.g., virus_A6
        meta = files.METADATA,
        seq = files.SEQUENCES


##############################
# Download from NBCI Virus with ingest snakefile
###############################
if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=files.SEQUENCES,
            metadata=files.METADATA
        threads: workflow.cores
        shell:
            """
            cd {input.dir}
            snakemake --cores {threads} all
            cd ../
            """

##############################
# AUGUR CURATE AND MERGE
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
# Merge with other metadata files you might have
###############################

rule curate:
    message:
        """
        Cleaning up metadata with augur curate
        """
    input:
        metadata = files.METADATA,  # Path to input metadata file
        meta_collab = files.meta_collab  # Data shared with us by collaborators
    params:
        strain_id_field=config["id_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        tmp = temp("temp/merged_meta.tsv"),  # Final output file for publications metadata
    output:
        meta="data/all_metadata.tsv"  # Final merged output file
    shell:
        """
        mkdir -p temp

        # Merge curated metadata
        augur merge --metadata metadata={input.metadata} meta_collab={input.meta_collab}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {params.tmp}

        # Normalize strings and format dates for metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {params.tmp} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.meta}
        echo "Curated metadata saved to {output.meta}"
        """

##############################
# BLAST
# blast fasta files for your specific proteins
# cut out your protein from fasta sequences
###############################

rule extract:
    input: 
        genbank_file = files.reference
    output: 
        extracted_fasta = "{seg}/config/reference.fasta",    
        extracted_genbank = "{seg}/config/reference.gbk",
        annotation = "{seg}/config/annotation.gff3"
    params:
        product_name = "{seg}",
        taxid = TAXID,
    shell:
        """
        python scripts/extract_gene_from_whole_genome.py \
        --genbank_file {input.genbank_file} \
        --output_fasta {output.extracted_fasta} \
        --product_name {params.product_name} \
        --output_genbank {output.extracted_genbank} \
        --taxid {params.taxid} \
        --output_gff {output.annotation}
        """

rule blast:
    input: 
        blast_db_file = rules.extract.output.extracted_fasta,  
        seqs_to_blast = rules.fetch.output.sequences
    output:
        blast_out = "temp/{seg}/blast_out.csv"
    params:
        blast_db =  "temp/{seg}/blast_database"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} \
        -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out, # output blast (for your protein)
        input_seqs = rules.fetch.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        range = "{seg}",  # Determines which protein (or whole genome) is processed - names must match seg in wildcard_constraints
        min_length = lambda wildcards: {"protein_xy": 600, "protein_yz": 2000, "genome": 6400}[wildcards.seg],  # Min length
        max_length = lambda wildcards: {"protein_xy": 900, "protein_yz": 2600, "genome": 8000}[wildcards.seg]  # Max length
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range} \
            --min_length {params.min_length} \
            --max_length {params.max_length}
        """

##############################
# Indexing sequences and filter them.
###############################

rule index_sequences:
    message:
        """
        Segment: {wildcards.seg}
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
        Segment: {wildcards.seg}

        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.curate.output.meta,
        exclude = files.exclude,
        include = files.include
    output:
        sequences = "{seg}/results/filtered.fasta",
        reason ="{seg}/results/reasons.tsv",
    log:
        "logs/filter.{seg}.log"
    params:
        group_by = "country year",
        sequences_per_group = 1000, # add a limit per group: correct for sampling bias
        strain_id_field= config["id_field"],
        min_date = 1950  # add a reasonable min date
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --include {input.include} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --output-sequences {output.sequences}\
            --output-log {output.reason} \
            >> {log} 2>&1

        echo "Filtered sequences saved to {output.sequences}"
        """
# --exclude-where ... or other parameters can be added, see `augur filter --h` for more options
##############################
# Reference for alignment added to sub-folders
###############################
rule align: 
    message:
        """
        Segment: {wildcards.seg}
        Aligning sequences to {input.reference} using Nextclade run.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = rules.filter.output.sequences,
        reference = rules.extract.output.extracted_fasta,

    output:
        alignment = "{seg}/results/aligned.fasta",
        tsv = "{seg}/results/nextclade.tsv",    
    params:
        penalty_gap_extend = config["align"]["penalty_gap_extend"],
        penalty_gap_open = config["align"]["penalty_gap_open"],
        penalty_gap_open_in_frame = config["align"]["penalty_gap_open_in_frame"],
        penalty_gap_open_out_of_frame = config["align"]["penalty_gap_open_out_of_frame"],
        kmer_length = config["align"]["kmer_length"],
        kmer_distance = config["align"]["kmer_distance"],
        min_match_length = config["align"]["min_match_length"],
        allowed_mismatches = config["align"]["allowed_mismatches"],
        min_length = config["align"]["min_length"]
    threads: workflow.cores
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations "{wildcards.seg}/results/translations/cds_{{cds}}.translation.fasta" \
        --output-fasta {output.alignment}
        """


##############################
# Tree building
###############################
rule tree:
    message:
        """
        Segment: {wildcards.seg}
        Creating a maximum likelihood tree with IQ-TREE
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
            --nthreads {threads} \
            --output {output.tree}
        """

##############################
# Refine to a timeline
###############################

rule refine:
    message:
        """
        Segment: {wildcards.seg}
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.curate.output.meta
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 3, # set to 6 if you want more control over outliers
        strain_id_field = config["id_field"],
        clock_rate = 0.004, # remove for estimation by augur; check literature
        clock_std_dev = 0.0015
    log:
        "logs/refine.{seg}.log"
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
            --stochastic-resolve \
            --clock-rate {params.clock_rate}\
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            2>&1 | (grep -i "pruning leaf" || cat > /dev/null) > {log}

        echo "Refined tree saved to {output.tree}"
        """
        #            

##############################
# Ancestral sequences and amino acids
###############################

rule ancestral:
    message: 
        """
        Segment: {wildcards.seg}

        Reconstructing ancestral sequences and mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
        annotation = rules.extract.output.extracted_genbank,
    output:
        node_data = "{seg}/results/muts.json",
    params:
        inference = "joint",
        genes = lambda wildcards: (
            CODING_GENES if wildcards.seg == "genome"
            else "protein_xy" if wildcards.seg == "protein_xy"
            else "protein_yz" if wildcards.seg == "protein_yz"
            else []
        ),        
        translation_template= r"{seg}/results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"{seg}/results/translations/cds_%GENE.ancestral.fasta",
        root = "{seg}/results/ancestral_sequences.fasta",
    log:
        "logs/ancestral.{seg}.log"
    shell:
        """
            augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template} \
            --output-sequences {params.root} \
            --inference {params.inference} \
            > {log} 2>&1
        """
        # --keep-ambiguous\ #do not infer nucleotides at ambiguous (N) sites on tip sequences (leave as N).
        # --root-sequence {input.annotation} \  -> assigns mutations to the root relative to the reference, not wanted here


rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata =rules.curate.output.meta
    output:
        node_data = "{seg}/results/traits.json" 
    params:
        traits = "country",
        strain_id_field= config["id_field"]
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
        muts = rules.ancestral.output.node_data,
        clades = files.clades # TODO: assign mutations to specific clades
    output:
        clade_data = "{seg}/results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.muts} \
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
        metadata = rules.curate.output.meta,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        muts = rules.ancestral.output.node_data,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= config["id_field"],
    output:
        auspice_json = "auspice/YOUR_VIRUS_{seg}.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

##############################

rule clean:
    message: 
        """
        Removing previous results and temporary files in the following directories:
        {params.targets}
        """
    params:
        targets = [
            "*/results/*",
            "auspice/*.json",
            "ingest/data/*.*",
            "temp/*",
            "logs/*",
            "ingest/benchmarks",
            "ingest/logs",
            "ingest/results",
            "benchmark/*"
            "data/fetch/*",
            "ingest/data/nextclade_data"
        ]
    shell: 
        """
        for dir in {params.targets}; do
            rm -rf "$dir" 2>/dev/null || true
        done
        """