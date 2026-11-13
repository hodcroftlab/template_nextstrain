# Nextstrain Template

This repository provides a reproducible workflow for building Nextstrain phylogenetic analyses of your virus. You can run analyses for specific proteins (e.g. VP1) and/or the whole genome.

For questions about Nextstrain or installation, refer to the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Checklist: What You Must Edit](#checklist-what-you-must-edit)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Repository Organization](#repository-organization)
- [Setup Instructions](#setup-instructions)
  - [Update Snakefile Parameters](#update-snakefile-parameters)
  - [Generate Reference Files](#generate-reference-files)
  - [Update Configuration Files](#update-configuration-files)
  - [Prepare Input Data](#prepare-input-data)
- [Adding a New Protein Segment](#adding-a-new-protein-segment)
- [Running Analyses](#running-analyses)
- [Visualizing Results](#visualizing-results)
- [Ingest Workflow](#ingest-workflow)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)

## Prerequisites

Ensure you have the following installed:

- Micromamba or Conda
- Python ≥ 3.8
- Snakemake ≥ 8 (tested with 9.x)
- Augur, Auspice, Nextclade CLI, BLAST+ (installed below)

The commands below assume Linux (Ubuntu/WSL).

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
      augur auspice "nextclade>=3" blast \
      snakemake git ncbi-datasets-cli unzip

    micromamba activate nextstrain
    ```

3. **Install additional dependencies:**

    ```bash
    micromamba install -c conda-forge -c bioconda csvtk seqkit tsv-utils ipdb entrez-direct
    micromamba install -c conda-forge fuzzywuzzy python-dotenv ipykernel
    ```

## Checklist: What You Must Edit

Before your first run, go through this list. Settings in the `snakefile` are referred to by name; use search (Ctrl+F) to find them.

| Where | What | Example |
| ----- | ---- | ------- |
| `snakefile` | `TAXID` | `"42789"` |
| `snakefile` | `YOUR_VIRUS` in output file names (replace all, Ctrl+H) | `ev_d68` |
| `snakefile` | `wildcard_constraints` → `seg` and the `segments` list | `"vp1\|genome"` and `['vp1', 'genome']` |
| `snakefile` | `min_length` / `max_length` per segment in rule `blast_sort` | `{"vp1": 600, "genome": 6400}` |
| `snakefile` | `CODING_GENES` (CDS names of the whole genome) | `["VP4", "VP2", ..., "3D"]` |
| `snakefile` | `genes` per segment in rule `ancestral` | `"VP1"` for `vp1` |
| `snakefile` | `clock_rate` / `clock_std_dev` in rule `refine` (check the literature) | `0.004` / `0.0015` |
| `snakefile` | `min_date`, `group_by`, `sequences_per_group` in rule `filter` | `1950` |
| `ingest/config/config.yaml` | `entrez_search_term`, `ncbi_taxon_id` | `"Enterovirus D68"`, `42789` |
| `<segment>/config/auspice_config.json` | Title, maintainers, colorings | |
| `<segment>/config/clades_genome.tsv` | Clade-defining mutations | |

## Quick Start

1. Complete the [checklist](#checklist-what-you-must-edit) above.
2. **Generate reference files** (see [Generate Reference Files](#generate-reference-files)).
3. **Run the workflow** (this also downloads the sequences via the ingest workflow):

    ```bash
    snakemake --cores 9 all
    ```

4. **View results:**

    ```bash
    auspice view --datasetDir auspice
    ```

## Pipeline Overview

The `snakefile` runs these steps for each segment:

1. **fetch** — runs the `ingest/` workflow to download sequences and metadata from NCBI into `data/fetch/`
2. **curate** — merges NCBI metadata with `data/meta_collab.tsv` and standardizes dates (`augur merge`, `augur curate`)
3. **extract** — cuts the segment (protein or whole genome) out of the GenBank reference
4. **blast** / **blast_sort** — finds the segment in each sequence with BLAST and keeps sequences within the length range
5. **index_sequences** / **filter** — subsamples and excludes sequences (`augur filter`)
6. **align** — aligns to the segment reference (`nextclade run`)
7. **tree** / **refine** — builds a maximum-likelihood tree and a time tree (`augur tree`, `augur refine`)
8. **ancestral** / **traits** / **clades** — reconstructs mutations, ancestral traits (country) and clades
9. **export** — writes the Auspice JSON to `auspice/`

## Repository Organization

- **`config/`** — Shared configuration files:
  - `config.yaml` — Main analysis parameters (ID column, date formats, alignment settings)
  - `colors.tsv` — Color scheme for tree nodes and traits
  - `geo_regions.tsv` — Geographic region definitions
  - `lat_longs.tsv` — Geographic coordinates
  - `exclude.txt` — Strain IDs to exclude from the analysis (one per line)
  - `include.txt` — Strain IDs to always keep, regardless of subsampling (one per line)
  - `reference.gbk`, `reference.fasta`, `genome_annotation.gff3` — Whole-genome reference files (you generate these, see [Generate Reference Files](#generate-reference-files))

- **`protein_xy/config/` and `genome/config/`** — Segment-specific configs:
  - `auspice_config.json` — Display settings for Auspice (colors, filters, defaults)
  - `clades_genome.tsv` — Clade/subgenotype definitions (nucleotide and amino acid mutations)
  - `reference.fasta`, `reference.gbk` — Segment reference (generated automatically by the workflow)
  - `annotation.gff3` — Segment annotation (generated automatically for protein segments; for `genome`, see [Generate Reference Files](#generate-reference-files))

- **`data/`** — Sequence and metadata files:
  - `fetch/sequences.fasta` — Sequences (from ingest or manual download)
  - `fetch/metadata.tsv` — Metadata (from ingest or manual download)
  - `meta_collab.tsv` — Optional collaborator metadata to merge
  - `all_metadata.tsv` — Merged and curated metadata (generated by the workflow)

- **`ingest/`** — Data download workflow (see [Ingest Workflow](#ingest-workflow))

- **`scripts/`** — Custom Python scripts:
  - `extract_gene_from_whole_genome.py` — Extracts a protein (or the whole genome) from the GenBank reference
  - `blast_sort.py` — Cuts the BLAST hit out of each sequence and filters by length

- **`snakefile`** — The entire computational pipeline, managed using Snakemake. Snakemake documentation can be found [here](https://snakemake.readthedocs.io/en/stable/).

- **`<segment>/results/`** — Analysis outputs (alignments, trees, node data JSONs)

- **`auspice/`** — Final Auspice JSON files for visualization

## Setup Instructions

### Update Snakefile Parameters

Edit the `snakefile` (search for each name):

- **`TAXID`** — NCBI Taxonomy ID for your virus.
- **`YOUR_VIRUS`** — Replace all occurrences (Ctrl+H) with your virus name, e.g. `ev_d68`. Output files are named `auspice/YOUR_VIRUS_<segment>.json`.
- **`wildcard_constraints`** — List your segments separated by `|`, e.g. `seg="vp1|genome"`.
- **`segments`** — The same segment names as a Python list, e.g. `['vp1', 'genome']`. These must match `wildcard_constraints` exactly.
- **`min_length` / `max_length`** in rule `blast_sort` — Expected length range (in nucleotides) of each segment. Sequences outside this range are dropped.
- **`CODING_GENES`** — Names of all CDS in your whole-genome reference. Used by `augur ancestral` for the `genome` segment.
- **`FETCH_SEQUENCES`** — `True` downloads data via ingest; set to `False` if you provide the data yourself (see [Prepare Input Data](#prepare-input-data)).

The segment name `genome` is reserved for the whole-genome analysis: the scripts treat it differently from protein segments.

### Generate Reference Files

Download your reference sequence from GenBank and generate FASTA, GenBank and annotation files, e.g. for EV-D68 (accession `AY426531`):

```bash
python3 ingest/bin/generate_from_genbank.py --reference "AY426531" --output-dir config/
```

When prompted for CDS annotation selection:

- Enter `[0]` for the first option
- Enter `[product]` to use product names, or leave blank for manual selection
- Enter `[2]` for the final selection

**Generated files:**

- `config/reference.gbk` — GenBank reference (input for the workflow)
- `config/reference.fasta` — Whole-genome FASTA
- `config/genome_annotation.gff3` — Whole-genome annotation

The whole-genome alignment needs an annotation file in the segment folder. Copy it there:

```bash
cp config/genome_annotation.gff3 genome/config/annotation.gff3
```

Segment-specific references for proteins (e.g. `protein_xy/config/reference.fasta`) are generated automatically during the workflow. The segment name must match the CDS `product` (or `gene`) name in `config/reference.gbk`.

### Update Configuration Files

- **`config/config.yaml`** — Set the ID column (`id_field`, default `accession`), the date formats found in your metadata, and adjust alignment parameters if needed (gap penalties, k-mer settings).

- **`<segment>/config/auspice_config.json`** — Customize display settings:
  - Title, maintainers, data provenance
  - Colorings (e.g., by clade, country, date)
  - Geographic resolutions
  - Default visualization options

- **`<segment>/config/clades_genome.tsv`** — Define clades using mutations (see [Nextstrain clade documentation](https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html)).

- **`config/exclude.txt` / `config/include.txt`** — Strain IDs to drop or always keep.

### Prepare Input Data

You can obtain sequences and metadata in two ways:

- **Automatic (default):** With `FETCH_SEQUENCES = True`, `snakemake all` runs the ingest workflow for you. You only need to configure it (see [Ingest Workflow](#ingest-workflow)).
- **Manual:** Download from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) (search for your virus or its Taxid), set `FETCH_SEQUENCES = False` in the `snakefile`, and save the files as:
  - `data/fetch/sequences.fasta`
  - `data/fetch/metadata.tsv`

Requirements for the metadata:

- An `accession` column (or the column set as `id_field` in `config/config.yaml`) matching the FASTA headers.
- A `date` column. Dates are converted to ISO 8601 (YYYY-MM-DD); all formats present must be listed under `expected_date_formats` in `config/config.yaml`. When saving from Excel, check that dates were not reformatted.

**Collaborator metadata (`data/meta_collab.tsv`)** is merged with the NCBI metadata. It must contain the same ID column (`accession`), otherwise `augur merge` fails. If you have no extra metadata, keep a file with only the header line.

## Adding a New Protein Segment

To analyse another protein, e.g. `vp1`:

1. **Create the config folder** by copying an existing one:

    ```bash
    cp -r protein_xy vp1
    ```

    Then update `vp1/config/auspice_config.json` and `vp1/config/clades_genome.tsv`.

2. **Use the same name everywhere** in the `snakefile`:
   - `wildcard_constraints` → `seg="vp1|genome"`
   - `segments` → `['vp1', 'genome']`
   - `min_length` / `max_length` in rule `blast_sort` → add `"vp1": ...`
   - `genes` in rule `ancestral` → add the CDS name for `vp1`

3. **Check the reference:** `vp1` must match the CDS `product` or `gene` name in `config/reference.gbk` (case-insensitive), otherwise the extraction fails.

Remove segments you do not use (e.g. the `protein_yz` example) from both `wildcard_constraints` and `segments`.

## Running Analyses

Activate the environment and run the workflow:

```bash
micromamba activate nextstrain
```

**Build all segments listed in `segments`:**

```bash
snakemake --cores 9 all
```

**Build a specific segment:**

```bash
snakemake auspice/YOUR_VIRUS_protein_xy.json --cores 9
snakemake auspice/YOUR_VIRUS_genome.json --cores 9
```

**Check what would run, without running it (dry run):**

```bash
snakemake -n all
```

**Start over:**

```bash
snakemake clean
```

Warning: `clean` deletes all results, **including the final Auspice JSONs in `auspice/`, the downloaded data in `data/fetch/`** and the ingest outputs. Back up anything you want to keep.

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

The `ingest/` subdirectory downloads and curates sequences from NCBI. It runs automatically as the `fetch` step of the main workflow, so you usually do not need to run it yourself.

**Configuration:**

- Edit `ingest/config/config.yaml` to set:
  - `entrez_search_term` — search query for your virus, e.g. `"Enterovirus D68"`
  - `ncbi_taxon_id` — NCBI taxonomy ID
  - `ncbi_datasets_fields` — metadata fields to retrieve

**Run ingest on its own (optional):**

```bash
cd ingest
snakemake --cores 9 all
cd ../
```

This writes:

- `data/fetch/sequences.fasta`
- `data/fetch/metadata.tsv`

For detailed ingest instructions, see [`ingest/README.md`](ingest/README.md).

The ingest pipeline is based on the Nextstrain [RSV ingest workflow](https://github.com/nextstrain/rsv.git).

## Troubleshooting

**Snakemake fails before running anything:**

- Run `snakemake -n all` to see the error.
- Check that every name in `segments` also appears in `wildcard_constraints` and has a `<segment>/config/` folder.

**Reference extraction fails ("Protein ... not found"):**

- Ensure `config/reference.gbk` exists and has complete CDS annotations with `product` or `gene` qualifiers.
- Check that the segment name matches one of those names.
- Use the interactive prompts in `generate_from_genbank.py` to select the correct features.

**Workflow fails at the alignment step:**

- For `genome`, check that `genome/config/annotation.gff3` exists (see [Generate Reference Files](#generate-reference-files)).
- Check that `min_length` / `max_length` in rule `blast_sort` fit your segment; too strict values leave no sequences.
- Verify alignment parameters in `config/config.yaml` are appropriate for your sequence diversity.

**Auspice doesn't load metadata correctly:**

- Confirm `data/fetch/metadata.tsv` has an `accession` column matching your FASTA sequence headers.
- Verify that all date formats in your metadata are listed in `config/config.yaml`.

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
