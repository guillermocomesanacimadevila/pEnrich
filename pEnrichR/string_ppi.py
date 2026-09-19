import stringdb
import polars as pl


class StringPPI:

    """
    Building a PPI network from the StringDB API 
    """

    def __init__(self):
        pass


    def query_string_enrichment(self, gene_list: list[str], add_nodes: int = 0) -> pl.DataFrame:
        gene_list = [i.upper() for i in gene_list]
        string_ids = stringdb.get_string_ids(gene_list)
        enrichment_df = stringdb.get_enrichment(string_ids.queryItem, add_nodes=add_nodes) ######### maybe I'll remove add_nodes!! 
        return pl.from_pandas(enrichment_df)


    def query_ppi(self, gene_list: list[str], add_nodes: int = 0) -> pl.DataFrame:
        gene_list = [i.upper() for i in gene_list]
        string_ids = stringdb.get_string_ids(gene_list)
        network_df = stringdb.get_network(string_ids.queryItem, add_nodes=add_nodes)
        network_df = pl.from_pandas(network_df)
        network_df = network_df.with_columns(
            pl.col("preferredName_A").is_in(set(gene_list)).alias("is_target_A"),
            pl.col("preferredName_B").is_in(set(gene_list)).alias("is_target_B"),
        )
        return network_df