# Nextstrain Template

This repository provides a reproducible workflow for building Nextstrain phylogenetic analyses of your virus. You can perform analyses for specific proteins or a full genome run.

For questions about Nextstrain or installation, refer to the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Repository Organization](#repository-organization)
- [Setup Instructions](#setup-instructions)
  - [Update Snakefile Parameters](#update-snakefile-parameters)
  - [Generate Reference Files](#generate-reference-files)
  - [Update Configuration Files](#update-configuration-files)
  - [Prepare Input Data](#prepare-input-data)
- [Running Analyses](#running-analyses)
- [Visualizing Results](#visualizing-results)
- [Ingest Workflow](#ingest-workflow)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)

## Prerequisites

Ensure you have the following installed:

- Python ≥ 3.8
- Micromamba or Conda
- Snakemake ≥ 7
- Nextstrain CLI

## Installation

1. **Clone the repository:**

    ```bash
    git clone git@github.com:hodcroftlab/template_nextstrain.git
    cd template_nextstrain
    ```

2. **Create and activate the Nextstrain environment:**

    ```bash
    micromamba create -n nextstrain \
      --override-channels --strict-channel-priority \
      -c conda-forge -c bioconda --yes \
      augur auspice nextclade \
      snakemake=7 git ncbi-datasets-cli

    micromamba activate nextstrain
    ```

3. **Install additional dependencies:**

    ```bash
    sudo apt-get update
    sudo apt-get install -y unzip

    micromamba install -c conda-forge -c bioconda csvtk seqkit tsv-utils ipdb entrez-direct
    micromamba install -c conda-forge fuzzywuzzy python-dotenv ipykernel
    ```

## Quick Start

1. **Set TAXID** in the `snakefile` (line 11)
2. **Generate reference files** (see [Setup Instructions](#setup-instructions))
3. **Run the workflow:**

    ```bash
    snakemake --cores 9 all
    ```

4. **View results:**

    ```bash
    auspice view --datasetDir auspice
    ```

## Repository Organization

This repository includes the following directories and files:

- **`config/`** — Shared configuration files:
  - `config.yaml` — Main analysis parameters (alignment settings, date formats)
  - `colors.tsv` — Color scheme for tree nodes and traits
  - `geo_regions.tsv` — Geographic region definitions
  - `lat_longs.tsv` — Geographic coordinates
  - `dropped_strains.txt` — Strain accessions to exclude from analysis
  - `reference_sequence.gb` — GenBank reference file (whole genome)

- **`protein_xy/config/` and `genome/config/`** — Segment-specific configs:
  - `auspice_config.json` — Display settings for Auspice (colors, filters, defaults)
  - `clades_genome.tsv` — Clade/subgenotype definitions (nucleotide and amino acid mutations)
  - `annotation.gff3` — Genome annotation (generated automatically)
  - `reference.fasta` — Segment reference sequence (generated automatically)

- **`data/`** — Sequence and metadata files:
  - `sequences.fasta` — Query sequences (from ingest or manual)
  - `metadata.tsv` — Sequence metadata (from ingest or manual)
  - `meta_collab.tsv` — Optional collaborator metadata to merge

- **`ingest/`** — Data ingestion workflow (see [Ingest Workflow](#ingest-workflow))

- **`scripts/`** — Custom Python scripts:
  - `extract_gene_from_whole_genome.py` — Extract proteins from GenBank reference
  - `blast_sort.py` — Sort and filter sequences by length

- **`snakefile`** — The entire computational pipeline, managed using Snakemake. Snakemake documentation can be found [here](https://snakemake.readthedocs.io/en/stable/).

- **`protein_xy/results/` and `genome/results/`** — Analysis outputs (alignments, trees, JSON data)

- **`auspice/`** — Final Auspice JSON files for visualization

## Setup Instructions

### Update Snakefile Parameters

Edit the `snakefile` to set your virus details:

- **Line 11:** Set `TAXID` to the NCBI Taxonomy ID for your virus
- **Line 30:** Set `wildcard_constraints` to match your analysis (e.g., `"vp1|whole_genome"`)
- **Line 34:** Update `segments` list to match your analysis (e.g., `['vp1', 'whole_genome']` - should match the above)
- **Line 37:** Define the minimum and maximum lengths of genes and the genome to include
- **Line 67:** Replace `<your_virus>` with your virus name in the output file naming

### Generate Reference Files

Extract your reference sequence (GenBank format) into FASTA and annotation files:

```bash
python3 ingest/bin/generate_from_genbank.py --reference "<accession>" --output-dir config/
```

When prompted for CDS annotation selection:

- Enter `[0]` for the first option
- Enter `[product]` to use product names, or leave blank for manual selection
- Enter `[2]` for the final selection

**Generated files:**

- `config/reference_sequence.gb` — GenBank reference
- `config/reference.fasta` — Whole genome FASTA (used by other rules)
- Segment-specific files are generated automatically during the workflow (e.g., `protein_xy/config/reference.fasta`)

### Update Configuration Files

- **`config/config.yaml`** — Adjust alignment parameters if needed (gap penalties, k-mer settings)

- **`protein_xy/config/auspice_config.json` and `genome/config/auspice_config.json`** — Customize display settings:
  - Title, maintainers, data provenance
  - Colorings (e.g., by clade, country, date)
  - Geographic resolutions
  - Default visualization options

- **`protein_xy/config/clades_genome.tsv` and `genome/config/clades_genome.tsv`** — Define clades using mutations (see [Nextstrain clade documentation](https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html))

### Prepare Input Data

You can obtain sequences and metadata via:

- **Automatic:** Run the ingest workflow (see [Ingest Workflow](#ingest-workflow))
- **Manual:** Download from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and save to:
  - `data/sequences.fasta`
  - `data/metadata.tsv`

Ensure metadata includes an `accession` column (or your configured ID field) and a `date` column in ISO format (YYYY-MM-DD).

## Running Analyses

Activate the environment and run the workflow:

```bash
micromamba activate nextstrain
```

**Build all segments (protein_xy + genome):**

```bash
snakemake --cores 9 all
```

**Build specific segment:**

```bash
snakemake auspice/<your_virus>_protein_xy.json --cores 9
snakemake auspice/<your_virus>_whole-genome.json --cores 9
```

**Clean intermediate files (keep final outputs):**

```bash
snakemake clean
```

## Visualizing Results

View your analyses locally with Auspice:

```bash
auspice view --datasetDir auspice
```

Open [http://localhost:4000](http://localhost:4000) in your browser.

**For simultaneous visualizations, set a different port:**

```bash
export PORT=4001
auspice view --datasetDir auspice
```

## Ingest Workflow

The `ingest/` subdirectory automates downloading and curating sequences from NCBI.

**Configuration:**

- Edit `ingest/config/config.yaml` to set:
  - `entrez_search_term` — search query for your virus, e.g. "Enterovirus D68"
  - `ncbi_taxon_id` — NCBI taxonomy ID
  - `ncbi_datasets_fields` — metadata fields to retrieve

**Run ingest:**

```bash
cd ingest
snakemake --cores 9 all
cd ../
```

This produces:

- `data/sequences.fasta`
- `data/metadata.tsv`

For detailed ingest instructions, see [`ingest/README.md`](ingest/README.md).

Sequences can be downloaded manually or automatically:

1. **Manual Download**: Visit [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), search for `<your_virus>` or Taxid `XXXXXX`, and download the sequences.
2. **Automated Download**: The `ingest` functionality handles automatic downloading via the workflow above.

The ingest pipeline is based on the Nextstrain [RSV ingest workflow](https://github.com/nextstrain/rsv.git).

## Troubleshooting

**Workflow fails at alignment step:**

- Ensure `config/reference_sequence.gb` exists and contains valid GenBank features
- Check that segment names (e.g., `protein_xy`, `genome`) match those in the `snakefile` line 34
- Verify alignment parameters in `config/config.yaml` are appropriate for your sequence diversity

**Auspice doesn't load metadata correctly:**

- Confirm `data/metadata.tsv` has an `accession` column matching your FASTA sequence headers
- Verify all dates are in ISO format (YYYY-MM-DD) using `config/config.yaml` date parameters

**Reference extraction fails:**

- Ensure GenBank file (`.gb`) has complete CDS annotations with `product` or `gene` qualifiers
- Use the interactive prompts in `generate_from_genbank.py` to select the correct features

For more help, see the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/) or the [Augur documentation](https://docs.nextstrain.org/projects/augur/en/stable/).

## Acknowledgments

- [Nextstrain](https://nextstrain.org/)
- [Auspice](https://auspice.us/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Augur](https://github.com/nextstrain/augur)
- [Biopython](https://biopython.org/)
- [NCBI](https://www.ncbi.nlm.nih.gov/)

## Contact

For questions or support, please open an [issue](https://github.com/hodcroftlab/template_nextstrain/issues) or contact the [hodcroftlab](https://eve-lab.org/).
