# pEnrichR

Pathway enrichment pipeline for pooled PPI (protein-protein interaction) networks.

Given a list of candidate gene targets, the pipeline:

1. **Queries [STRING](https://string-db.org/)** for the PPI network among the targets and filters edges by confidence score.
2. **Runs Markov Clustering (MCL)** on the filtered network to group genes into functional clusters.
3. **Runs [EnrichR](https://maayanlab.cloud/Enrichr/)** (GO Biological Process + KEGG) on each cluster.

Results are written per phenotype under `results/<pheno_id>/`.

## Installation

Requires Python >= 3.12.

```bash
python3.12 -m venv .venv
./.venv/bin/pip install -e .
```

> If your shell auto-activates a conda base environment, make sure `python3.12` above actually resolves to a 3.12 interpreter (`which python3.12`),  otherwise create the venv from an explicit interpreter path instead, e.g. `~/.pyenv/versions/3.12.8/bin/python3 -m venv .venv`.

## Usage

Once installed, run it as a console command:

```bash
.venv/bin/activate
```

```bash
pEnrichR --targets_file assets/targets.txt --pheno_id AD
```

or as a module/script without installing:

```bash
./.venv/bin/python main.py --targets_file assets/targets.txt --pheno_id AD
```

### Arguments

| Argument         | Required | Default   | Description                                      |
|------------------|----------|-----------|---------------------------------------------------|
| `--targets_file` | yes      | —         | Path to a text file of gene symbols, one per line |
| `--pheno_id`     | yes      | —         | Phenotype/run identifier, used to name outputs    |
| `--string_score` | no       | `0.4`     | Minimum STRING confidence score to keep an edge   |
| `--out_dir`      | no       | `results` | Output directory                                  |

### Output

```
results/<pheno_id>/
├── string_ppi/
│   └── <pheno_id>_string_ppi_filtered.tsv
├── mcl/
│   ├── <pheno_id>_mcl_clusters.png
│   └── <pheno_id>_mcl_clusters.tsv
└── enrichr/
    ├── <pheno_id>_cluster0_GO.tsv
    ├── <pheno_id>_cluster0_KEGG.tsv
    └── ...
```