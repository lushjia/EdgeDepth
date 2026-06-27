"""
Step 1 of redundancy filtering via the bubble-puncture method.

Given all edges in the pangenome graph and the GFA file for one chromosome,
this script does two things:

  1. Identifies bridge edges (edges whose removal would disconnect the graph).
     Bridge edges that are NOT on the GRCh38 reference walk are written to the
     bridge-edge output file.  These are kept as-is in the downstream filter.

  2. Decomposes the remaining graph (after removing all bridge edges) into
     biconnected components (at articualtion points).  Each large component (>= --min_nodes nodes) is
     saved as a separate GML file for independent processing in Step 2.
     All remaining small components are saved together as a single residual GML.

Usage:
    python filter_redundancy_bubble_puncture.step1_subgraph.py \\
        --chr          chr21 \\
        --edges        hprc-v2.0-mc-grch38.edges.txt \\
        --gfa          chr21.gfa \\
        --bridge_out   chr21.original_bridge_edge.txt \\ # before decompose graph, edges that not in GRCh38 edges, and independent
        --subgraph_dir chr21_subgraphs/

Output files written to --subgraph_dir:
    subgraph_1.graph.gml, subgraph_2.graph.gml, ...   (large biconnected components)
    subgraph_remaining.graph.gml                        (all small components merged)

"""

import sys 
import re 
import argparse
import os

import networkx as nx

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Identify bridge edges and decompose graph into biconnected components."
)
parser.add_argument(
    "--chr", required=True,
    help="chromosome to analyze"
)
parser.add_argument(
    "--edges", required=True,
    help="Edge list file (chr + edge_id columns, one edge per line, with header)."
)
parser.add_argument(
    "--gfa", required=True,
    help="GFA file for the chromosome; GRCh38 reference walk is read from W lines."
)
parser.add_argument(
    "--bridge_out", required=True,
    help="Output file for bridge edges that are not on the GRCh38 reference walk."
)
parser.add_argument(
    "--subgraph_dir", required=True,
    help="Output directory for biconnected-component GML files."
)
parser.add_argument(
    "--min_nodes", type=int, default=100,
    help="Biconnected components with fewer nodes than this are merged into the "
         "residual graph instead of being saved as individual files (default: 100)."
)
args = parser.parse_args()

os.makedirs(args.subgraph_dir, exist_ok=True)


# ---------------------------------------------------------------------------
# Build graph from edge file
# ---------------------------------------------------------------------------

# create graph from edge id 
G = nx.Graph()
node_edge_dict = {} # node1, node2: edge

# construct graph from all edges read depth raw count file 
with open(args.edges, 'r') as edge_raw_file: 
    next(edge_raw_file)  # skip header
    # chr edge_id
    for line in edge_raw_file:
        line = line.strip().split('\t')
        edge = line[1]
        if line[0] == args.chr:
            nodes = re.split(">|<", edge)[1:]
            # smaller node first 
            if nodes[0] > nodes[1]:
                nodes[0], nodes[1] = nodes[1], nodes[0]
            G.add_edge(nodes[0], nodes[1])
            node_edge_dict[tuple([nodes[0], nodes[1]])] = edge 

# ---------------------------------------------------------------------------
# Read GRCh38 reference edges from GFA walk lines
# ---------------------------------------------------------------------------
grch38_edge_set = set()  # to store GRCh38 edges # 102643231-102643232 ... 

with open(args.gfa, 'r') as gfa_file:
    for line in gfa_file:
        if not line.startswith('W'):
            continue
        parts = line.strip().split('\t')
        if parts[1] != "GRCh38":
            continue
        path = parts[6]
        # path is like: >102643231>102643232>102643233>102643234>102643235
        path_edges = re.split(">|<", path)[1:]
        # recording all edges in the path
        for i in range(len(path_edges) - 1):
            # smaller node first 
            node1 = path_edges[i]
            node2 = path_edges[i+1]
            if node1 > node2:
                node1, node2 = node2, node1
            edge = f"{node1}-{node2}"
            grch38_edge_set.add(edge)


# ---------------------------------------------------------------------------
# Find bridge edges and write non-reference ones to output
# ---------------------------------------------------------------------------
# get all brdige edges from origianl graph 
ori_bridge_edges = set(nx.bridges(G))  # set of tuples (node1, node2) # len 11886

# independent edge file output 

# ori bridge edges - independent expect reference edges 
with open(args.bridge_out, 'w') as original_bridge_edge_file:
    for edge in ori_bridge_edges:
        # pop edge from graph 
        G.remove_edge(*edge)
        # write into independent edge file if not in grch38_edge_set
        node1, node2 = edge
        # smaller node first
        if node1 > node2:
            node1, node2 = node2, node1
        edge_mod = f"{node1}-{node2}"
        if edge_mod not in grch38_edge_set:
            # write to file 
            edge_ori = node_edge_dict[(node1, node2)]
            original_bridge_edge_file.write(f"{edge_ori}\n")  # write original edge format >node1>node2
            # original_bridge_edge_file.flush()


# ---------------------------------------------------------------------------
# Decompose remaining graph into biconnected components
# ---------------------------------------------------------------------------
# copy graph to avoid modifying the original graph
G_copy = G.copy()  # create a copy of the graph to avoid modifying the original graph

n = 1
for bicomponent in nx.biconnected_components(G_copy):
    subgraph = G.subgraph(bicomponent)
    subgraph = nx.Graph(subgraph)  # unfreeze the graph
    # if nodes number <50, then keep in graph and save as a whole later 
    if len(subgraph.nodes) < args.min_nodes:
        continue  # skip small subgraphs
    # save into a gml file 
    out_path = os.path.join(args.subgraph_dir, f"subgraph_{n}.graph.gml")
    nx.write_gml(subgraph, out_path)
    # nx.write_gml(subgraph, f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/bubble_puncture/{chr}/subgraph_{n}.graph.gml")
    ebunch = list(subgraph.edges())
    G.remove_edges_from(ebunch)
    n += 1
    # read it back 
    # G_loaded = nx.read_gml("graph.gml")

# Save residual (small components and any remaining edges)
# remove the isolated nodes
G.remove_nodes_from(list(nx.isolates(G)))
# save the remain graph as a whole
nx.write_gml(G, os.path.join(args.subgraph_dir, "subgraph_remaining.graph.gml"))



