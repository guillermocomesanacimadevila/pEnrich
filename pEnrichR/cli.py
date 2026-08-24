import argparse
import polars as pl
from .string_ppi import StringPPI
from .unsupervised_algorithms import MarkovClustering
from .enrichr import EnrichR
from .utils import (
    pathway_string_ppi_out,
    pathway_mcl_plot_out,
    pathway_mcl_clusters_out,
    pathway_enrichr_out,
    txt_to_list
)

def query_and_filter_ppi(targets: list[str], string_score: float, pheno_id: str, out_dir: str = "results"):

    """
    Query STRING PPI for the given targets and filter by score.
    """

    ppi_network = StringPPI().query_ppi(targets) # just a df for STRING results
    # add score gait into params/ for STRING ##################
    ppi_network = ppi_network.filter(pl.col("score").cast(pl.Float64) > string_score)
    for row in ppi_network.iter_rows(named=True):
        gene1 = row["preferredName_A"]
        gene2 = row["preferredName_B"]
        sc = row["score"]
        print(f"[TRACKING] STRING score between {gene1} and {gene2}: {sc}")

    out_file = pathway_string_ppi_out(pheno_id, out_dir)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    ppi_network.write_csv(str(out_file), separator="\t")
    print(f"[DONE] Saved filtered STRING PPI: {out_file}")
    return ppi_network


def build_ppi_network(targets_file: str, pheno_id: str, string_score: float, out_dir: str = "results"):

    """
    Grab candidate targets from a text file -> build + filter STRING PPI.
    """

    targets = txt_to_list(targets_file)

    return query_and_filter_ppi(
        targets=targets,
        string_score=string_score,
        pheno_id=pheno_id,
        out_dir=out_dir,
    )


def perform_mcl_clustering(string_ppi: pl.DataFrame, pheno_id: str, out_dir: str = "results"):

    """
    MCL clustering
    """

    mcl = MarkovClustering()

    matrix = mcl.df_to_matrix(
        gene_1_col=string_ppi["preferredName_A"],
        gene_2_col=string_ppi["preferredName_B"],
        score_col=string_ppi["score"]
    )

    result, index_clusters, gene_clusters, Q = mcl.markov_clustering(matrix=matrix)

    for cluster in gene_clusters:
        print(f"[TRACKING] MCL cluster ({len(cluster)} genes): {cluster}")

    plot_file = pathway_mcl_plot_out(pheno_id, out_dir)
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    mcl.draw_clusters(matrix, index_clusters, output_path=str(plot_file))
    print(f"[DONE] Saved MCL cluster plot: {plot_file}")

    clusters_file = pathway_mcl_clusters_out(pheno_id, out_dir)
    clusters_df = pl.DataFrame({
        "cluster_id": [i for i, cluster in enumerate(gene_clusters) for _ in cluster],
        "gene": [gene for cluster in gene_clusters for gene in cluster],
    })
    clusters_df.write_csv(str(clusters_file), separator="\t")
    print(f"[DONE] Saved MCL clusters ({len(gene_clusters)} clusters, modularity Q={Q}): {clusters_file}")
    return gene_clusters


def run_enrichr(cluster_list: list[list], pheno_id: str, out_dir: str = "results"):

    """
    For each cluster -> run EnrichR
    """

    enrichr = EnrichR()
    results = []
    for cluster_id, cluster in enumerate(cluster_list):
        go_df, kegg_df = enrichr.run_enrichr(gene_list=cluster)
        print(f"[TRACKING] EnrichR done for cluster ({len(cluster)} genes): {cluster}")
        go_file = pathway_enrichr_out(pheno_id, cluster_id, "GO", out_dir)
        go_file.parent.mkdir(parents=True, exist_ok=True)
        go_df.write_csv(str(go_file), separator="\t")
        kegg_file = pathway_enrichr_out(pheno_id, cluster_id, "KEGG", out_dir)
        kegg_df.write_csv(str(kegg_file), separator="\t")
        print(f"[DONE] Saved EnrichR results for cluster {cluster_id}: {go_file}, {kegg_file}")
        results.append((cluster, go_df, kegg_df))
    return results


def main():
    p = argparse.ArgumentParser(prog="pEnrichR")
    p.add_argument("--targets_file", required=True)
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--string_score", type=float, default=0.4)
    p.add_argument("--out_dir", default="results")
    args = p.parse_args()

    ppi_network = build_ppi_network(
        targets_file=args.targets_file,
        pheno_id=args.pheno_id,
        string_score=args.string_score,
        out_dir=args.out_dir,
    )

    gene_clusters = perform_mcl_clustering(
        ppi_network,
        pheno_id=args.pheno_id,
        out_dir=args.out_dir,
    )

    results = run_enrichr(
        gene_clusters,
        pheno_id=args.pheno_id,
        out_dir=args.out_dir,
    )

    for cluster, go_df, kegg_df in results:
        print(f"[TRACKING] cluster: {cluster}")
        print(go_df)
        print(kegg_df)


if __name__ == "__main__":
    main()
