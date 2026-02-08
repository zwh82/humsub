# Subspecies Quantification Workflow – `use_catalog`

This Snakemake pipeline quantifies microbial **subspecies** in metagenomic samples using `sourmash gather` and the HuMSub catalog.

It is part of the [HuMSub](https://zenodo.org/records/15862096) framework for high-resolution metagenomic analysis.

### python

```
mamba create -n gut sourmash=4.6.1 python=3.9 path pandas pyyaml -y

# download
python humsub/Workflows/use_catalog/workflow/python/main.py --task download --config configs/default_config_python.yaml

# process 
python humsub/Workflows/use_catalog/workflow/python/main.py --task process --config configs/default_config_python.yaml --output-dir output --threads 16

# combine
python humsub/Workflows/use_catalog/workflow/python/main.py --task combine --config configs/default_config_python.yaml --output-dir output --threads 16

# taxonomy
python humsub/Workflows/use_catalog/workflow/python/main.py --task taxonomy --config configs/default_config_python.yaml --output-dir output --threads 16
```

---

## Overview

The workflow:

1. (Optionally) generates a `samples.tsv` sample sheet from raw FASTQ files  
2. Optionally trims input reads (Check [link](https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html#how-should-i-prepare-my-data) for details)
3. Creates sourmash sketches (`.sig`) from each sample  
4. Performs `sourmash gather` against a prebuilt subspecies SBT  
5. Aggregates subspecies hits and maps them to taxonomy

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- `conda` or `mamba` for environment management
- FASTQ files or a pre-defined `samples.tsv`

---

## Setup

```bash
# Clone the full project
git clone https://github.com/trajkovski-lab/humsub.git
cd humsub/Workflows/use_catalog

# Copy and edit config
cp ../config/default_config.yaml config.yaml
# Edit paths inside config.yaml to match your setup
```

## Usage

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile config.yaml \
          --use-conda \
          --cores 4
```

## Inputs

You can either:

- Place your raw `.fastq.gz` files in a directory and set `generate_sample_table: true` in the config  
- Provide a custom `samples.tsv` file with sample metadata (must include `Reads_QC_R1`, `Reads_QC_R2`, or `Reads_QC_se` columns, please use [metagenome-atlas](https://github.com/metagenome-atlas/atlas) init command.)

These options are controlled through `config.yaml`.

---

## Outputs

- `output/subspecies_taxonomy.csv` – Final output: subspecies relative abundances mapped to taxonomy  
- `output/subspecies_relab.csv` – Intermediate file: raw relative abundances per sample  
- `output/sketch_samples/` – Sourmash signature (`.sig`) files  
- `output/gather/` – Per-sample gather CSVs (unaggregated)

---

## License

This workflow is distributed under the terms specified in the root `LICENSE.txt`.