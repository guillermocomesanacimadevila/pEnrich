from pathlib import Path
import os 

def pathway_string_ppi_out(pheno_id: str, out_dir: str | Path = "results") -> Path:
    return (Path(out_dir) / pheno_id / "string_ppi" / f"{pheno_id}_string_ppi_filtered.tsv")


def pathway_mcl_plot_out(pheno_id: str, out_dir: str | Path = "results") -> Path:
    return (Path(out_dir) / pheno_id / "mcl" / f"{pheno_id}_mcl_clusters.png")


def pathway_mcl_clusters_out(pheno_id: str, out_dir: str | Path = "results") -> Path:
    return (Path(out_dir) / pheno_id / "mcl" / f"{pheno_id}_mcl_clusters.tsv")


def pathway_mcl_sweep_out(pheno_id: str, out_dir: str | Path = "results") -> Path:
    return (Path(out_dir) / pheno_id / "mcl" / f"{pheno_id}_mcl_inflation_sweep.tsv")


def pathway_enrichr_out(pheno_id: str, cluster_id: int, library: str, out_dir: str | Path = "results") -> Path:
    return (Path(out_dir) / pheno_id / "enrichr" / f"{pheno_id}_cluster{cluster_id}_{library}.tsv")


def txt_to_list(txt_file):
    ls = []
    with open(txt_file, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                ls.append(line)
    return set(ls)

