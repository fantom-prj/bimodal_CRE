# Network basic statistics (node/edge composition across variants)

Regenerates node and edge counts across all 8 network variants (4 base + 4
lasso).

* `build_network_composition_data.R` — after applying the significance
  filter to each of the 8 candidate network variants (base vs elastic-net-
  regularised, with/without TOBIAS footprinting, with/without Hi-C
  contacts), counts distinct nodes (TFs, enhancers, promoters, by
  measurement type) and edges (by link type and measurement-type
  combination). Reads all 8 network RDS files via a `primary_data_folder`
  placeholder.

Writes `Fig5/data/network_node_counts.csv` and
`Fig5/data/network_edge_counts.csv` (both already display-labelled: network
as an ordered factor, metric as an ordered factor, plus a `colour` column for
direct ggplot use).

Feeds `Fig5/ext_fig5c_network_nodes.R` and `Fig5/ext_fig5c_network_edges.R`.
