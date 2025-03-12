# Nextstrain Template

This repository provides a comprehensive Nextstrain analysis of "your virus". You can choose to perform either a shorter run with specific proteins or a full genome run.

For those unfamiliar with Nextstrain or needing installation guidance, please refer to the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/).

## Table of Contents
- [Prerequisites](#prerequisites)
- [Nextstrain Environment](#nextstrain-environment)
  - [Installation](#installation)
  - [Enhancing the Analysis](#enhancing-the-analysis)
- [Repository Organization](#repository-organization)
  - [Configuration Files](#configuration-files)
- [Usage Examples](#usage-examples)
  - [Running a Build](#running-a-build)
  - [Visualizing the Build](#visualizing-the-build)
- [Ingest](#ingest)
  - [Sequences](#sequences)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)

## Prerequisites
Ensure you have the following installed:
- Python=3.8 or higher
- Micromamba or Conda
- Snakemake=7
- Nextstrain CLI

## Nextstrain Environment
Install the Nextstrain environment by following [these instructions](https://docs.nextstrain.org/en/latest/guides/install/local-installation.html).

### Installation
1. Clone the repository:
    ```bash
    git clone git@github.com:hodcroftlab/template_nextstrain.git
    cd template_nextstrain
    ```

2. Install the Nextstrain environment:
    ```bash
    micromamba create -n nextstrain \
      --override-channels --strict-channel-priority \
      -c conda-forge -c bioconda --yes \
      augur auspice nextclade \
      snakemake=7 git ncbi-datasets-cli

    micromamba activate nextstrain
    ```

3. Update/install additional dependencies:
    ```bash
    sudo apt-get update
    sudo apt-get install -y unzip

    micromamba install -c conda-forge -c bioconda csvtk seqkit tsv-utils ipdb entrez-direct
    micromamba install -c conda-forge fuzzywuzzy python-dotenv ipykernel
    ```

### Enhancing the Analysis
The data for this analysis is available from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). Instructions for downloading sequences are provided under [Sequences](#sequences).

## Repository Organization
This repository includes the following directories and files:

- `scripts`: Custom Python scripts called by the `snakefile`.
- `snakefile`: The entire computational pipeline, managed using Snakemake. Snakemake documentation can be found [here](https://snakemake.readthedocs.io/en/stable/).
- `ingest`: Contains Python scripts and the `snakefile` for automatic downloading of <your_virus> sequences and metadata.
- <`protein_xy`>: Sequences and configuration files for the specific **protein_xy run**.
- `whole_genome`: Sequences and configuration files for the **whole genome run**.

### Configuration Files
The `config`, `protein_xy/config`, and `whole_genome/config` directories contain necessary configuration files:
- `config.yaml`: Configuration file for setting parameters and options for the analysis
- `colors.tsv`: Color scheme
- `geo_regions.tsv`: Geographical locations
- `lat_longs.tsv`: Latitude data
- `dropped_strains.txt`: It will exclude these accessions during `augur filter`
- `clades_genome.tsv`: Manually Labeling Clades on a Nextstrain Tree (see documentation [here](https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html))
- `reference_sequence.gb`: Reference sequence (add manually)
- `auspice_config.json`: Auspice configuration file - has to be in all data folders!

The reference sequence used is [XYZ, accession number](https://www.ncbi.nlm.nih.gov/nucleotide/), sampled in 19XX.

## Usage Examples
### Running a Build
Activate the Nextstrain environment:
```bash
micromamba activate nextstrain
```

To perform a build, run:
```bash
snakemake --cores 9 all
```

For specific builds:
- protein_xy build:
```bash
snakemake auspice/<your_virus>_protein_xy.json --cores 9
```
- Whole genome build:
```bash
snakemake auspice/<your_virus>_whole-genome.json --cores 9
```

### Visualizing the Build
To visualize the build, use Auspice:
```bash
auspice view --datasetDir auspice
```
To run two visualizations simultaneously, you may need to set the port:
```bash
export PORT=4001
```

## Ingest
For more information on how to run the `ingest`, please refer to the [README](ingest/README.md) in the `ingest` folder.

### Sequences
Sequences can be downloaded manually or automatically.

1. **Manual Download**: Visit [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), search for `<your_virus>` or Taxid `XXXXXX`, and download the sequences.
2. **Automated Download**: The `ingest` functionality, included in the main `snakefile`, handles automatic downloading.

The ingest pipeline is based on the Nextstrain [RSV ingest workflow](https://github.com/nextstrain/rsv.git). Running the **ingest** pipeline produces `data/metadata.tsv` and `data/sequences.fasta`.


## Acknowledgments
- [Nextstrain](https://nextstrain.org/)
- [Auspice](https://auspice.us/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Biopython](https://biopython.org/)
- [Genbank](https://www.ncbi.nlm.nih.gov/genbank/)
- [NCBI](https://www.ncbi.nlm.nih.gov/)

## Contact
For questions or support, please contact [your-email@example.com].

Feel free to adjust the content according to your project's specifics.