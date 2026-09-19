import markov_clustering as mc
import networkx as nx
import polars as pl
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import re


class MarkovClustering:

    """
    MCL (network-based) clustering for STRING-derived PPIs
    """

    def __init__(self):
        self.network = nx.Graph()
        self.node_order = []


    def df_to_matrix(self, gene_1_col: list[str], gene_2_col: list[str], score_col: list[float]):

        """
        For each row in df (i.e STRING results) - asumming they're filtered by conf,
        grab each pair of proteins (i.e. network walk) and  create a df ->
        pl.DataFrame({
            "nameA": ["GENE1", "GENE2"],
            "nameB": ["GENEX", "GENEY"],
            "score": ["1", "2"]
        })
        """

        network_df = pl.DataFrame({
            "preferredName_A": gene_1_col,
            "preferredName_B": gene_2_col,
            "score": score_col
        })

        for row in network_df.iter_rows(named=True):
            self.network.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

        self.node_order = list(self.network.nodes())
        matrix = sp.csr_matrix(nx.to_scipy_sparse_array(self.network, nodelist=self.node_order, weight="weight", format="csr"))
        return matrix
    

    def markov_clustering(self, matrix, inflation: float = 2.0):

        """
        Run MCL and convert the resulting matrix indices back into gene names.
        """

        result = mc.run_mcl(matrix, inflation=inflation)
        index_clusters = mc.get_clusters(result)
        gene_clusters = [[self.node_order[index] for index in cluster] for cluster in index_clusters]
        Q = mc.modularity(matrix=result, clusters=index_clusters)
        print("inflation:", inflation, "modularity:", Q)
        return result, index_clusters, gene_clusters, Q


    def sweep_inflation(self, matrix, inflation_values: list[float] = [1.4, 2.0, 2.5, 3.0, 4.0]):

        """
        Run MCL across a range of inflation values, pick the one with highest modularity Q.
        """

        rows = []
        best = None
        for inflation in inflation_values:
            result, index_clusters, gene_clusters, Q = self.markov_clustering(matrix, inflation=inflation)

            rows.append({
                "inflation": inflation,
                "n_clusters": len(gene_clusters),
                "modularity_q": Q,
            })

            if best is None or Q > best[-1]:
                best = (inflation, result, index_clusters, gene_clusters, Q)

        sweep_df = pl.DataFrame(rows)
        best_inflation, result, index_clusters, gene_clusters, Q = best
        print(f"[DONE] Best inflation: {best_inflation} (modularity Q={Q})")
        return result, index_clusters, gene_clusters, Q, sweep_df


    def _packed_layout(self, seed: int = 42):

        """
        Lay out the main connected component on its own (with room to breathe), then
        pack every smaller disconnected component into a tight grid beside it - keeps
        small pairs from drifting into empty space the way a plain spring_layout does
        on a disconnected graph.
        """

        components = sorted(nx.connected_components(self.network), key=len, reverse=True)
        main_component, small_components = components[0], components[1:]

        main_subgraph = self.network.subgraph(main_component)
        k = 2.6 / np.sqrt(len(main_component))
        pos = nx.spring_layout(main_subgraph, seed=seed, k=k, iterations=800)

        xs = [x for x, y in pos.values()]
        ys = [y for x, y in pos.values()]

        n_cols = max(1, int(np.ceil(np.sqrt(len(small_components)))))
        cell = 0.45
        start_x = max(xs) + cell
        start_y = max(ys)

        for i, component in enumerate(small_components):
            row, col = divmod(i, n_cols)
            cx = start_x + col * cell * 2.0
            cy = start_y - row * cell * 2.0

            if len(component) == 1:
                pos[next(iter(component))] = (cx, cy)
                continue

            if len(component) == 2:
                node_a, node_b = component
                pos[node_a] = (cx, cy)
                pos[node_b] = (cx + cell * 0.7, cy)
                continue

            subgraph = self.network.subgraph(component)
            sub_pos = nx.spring_layout(subgraph, seed=seed, k=0.3 / np.sqrt(len(component)), iterations=200, scale=cell * 0.4)
            sub_xs = [x for x, y in sub_pos.values()]
            sub_ys = [y for x, y in sub_pos.values()]
            for node, (x, y) in sub_pos.items():
                pos[node] = (x - min(sub_xs) + cx, y - min(sub_ys) + cy)

        return pos


    @staticmethod
    def _short_label(label: str) -> str:
        stripped = re.sub(r"(?i)^(positive |negative )?regulation of ", "", label)
        return stripped


    def draw_clusters(self, index_clusters, cluster_labels: dict[int, str], targets: set[str] = None, output_path: str = "string_mcl_clusters.png"):

        """
        Draw the clustered STRING network - nodes sized by degree, colored by cluster pathway label,
        dashed outline for genes not in targets (i.e. added via --add_nodes).
        """

        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]

        pos = self._packed_layout()
        degrees = dict(self.network.degree())

        unique_labels = list(dict.fromkeys(cluster_labels.values()))
        cmap = plt.get_cmap("tab10")
        label_colors = {label: cmap(i % 10) for i, label in enumerate(l for l in unique_labels if l != "No significant pathway")}
        label_colors["No significant pathway"] = (0.82, 0.82, 0.82, 1.0)

        fig, ax = plt.subplots(figsize=(13, 9))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        nx.draw_networkx_edges(self.network, pos, edge_color="#c9c9c9", width=0.8, alpha=0.7, ax=ax)

        for cluster_id, gene_indices in enumerate(index_clusters):
            genes = [self.node_order[index] for index in gene_indices]
            label = cluster_labels.get(cluster_id, "No significant pathway")
            for linestyle in ("solid", "dashed"):
                sub_genes = [
                    gene for gene in genes
                    if (targets is None or gene in targets) == (linestyle == "solid")
                ]
                if not sub_genes:
                    continue
                sizes = [70 + degrees[gene] * 30 for gene in sub_genes]
                collection = nx.draw_networkx_nodes(
                    self.network, pos, nodelist=sub_genes, node_color=[label_colors[label]] * len(sub_genes),
                    node_size=sizes, edgecolors="#222222", linewidths=1.0, ax=ax
                )
                collection.set_linestyle(linestyle)

        for gene, (x, y) in pos.items():
            ax.text(x, y + 0.045, gene, fontsize=6.5, ha="center", va="bottom", color="#111111")

        legend_handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markeredgecolor="#222222",
                       markersize=8, label=self._short_label(label))
            for label, color in label_colors.items()
        ]
        legend = ax.legend(handles=legend_handles, title="Enriched pathway", loc="best",
                            frameon=True, fontsize=7.5, title_fontsize=8.5, fancybox=True, framealpha=0.95)
        legend.get_frame().set_edgecolor("#999999")
        legend.get_frame().set_linewidth(0.7)

        ax.axis("off")
        ax.set_aspect("equal")
        ax.margins(0.08)
        plt.savefig(output_path, dpi=400, bbox_inches="tight", facecolor="white")
        plt.close()