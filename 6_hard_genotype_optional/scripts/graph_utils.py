"""
Shared graph helpers for the 6_hard_genotype_optional scripts.

Only depends on the standard library (no networkx/scipy), so it can be
imported from any step's container.
"""

import gzip
import re


def open_maybe_gzip(path):
    """Open a plain-text or .gz file for reading in text mode."""
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


# ---------------------------------------------------------------------------
# Handle / path / edge string helpers
#   handle : '>123' or '<123' (node id with traversal orientation)
#   path   : '>1>3<4' (concatenated handles, as in GFA W-lines)
#   edge   : two concatenated handles, e.g. '>13>14'
# ---------------------------------------------------------------------------
HANDLE_PATTERN = re.compile(r"([<>])([^<>]+)")


def parse_path(path):
    """Parse '>1>3<4' into ['>1', '>3', '<4']."""
    path = path.replace(" ", "")
    return [direction + node for direction, node in HANDLE_PATTERN.findall(path)]


def reverse_handle(handle):
    """Flip a handle's orientation: '>5' <-> '<5'."""
    return ("<" if handle[0] == ">" else ">") + handle[1:]


def reverse_complement_path(path):
    """Reverse-complement a path string, e.g. '>5>6<7' -> '>7<6<5'."""
    return "".join(reverse_handle(h) for h in reversed(parse_path(path)))


def canonical_edge(left, right):
    """
    Canonical string for an edge, treating an edge and its reverse-complement
    traversal as the same edge (e.g. '>13>14' == '<14<13').
    The representation starting with '>' is preferred.
    """
    candidates = [(left, right), (reverse_handle(right), reverse_handle(left))]
    chosen = min(
        candidates,
        key=lambda e: (e[0][0] != ">", e[0][1:], e[1][0], e[1][1:]),
    )
    return chosen[0] + chosen[1]


def precompute_component_traversal(genotype_graph):
    """
    Build the sample-independent traversal setup for one component's
    (condensed) directed multigraph, so it is computed once per chromosome
    instead of once per sample.

    Returns:
        incoming        : node -> [(source, target, key), ...] edges into node
        outgoing        : node -> [(source, target, key), ...] edges out of node
        edge_depth_keys : (source, target, key) -> canonical edge string used
                          to look the edge up in the depth matrix
    """
    incoming = {}
    outgoing = {}
    edge_depth_keys = {}
    for source, target, key in genotype_graph.edges(keys=True):
        edge = (source, target, key)
        outgoing.setdefault(source, []).append(edge)
        incoming.setdefault(target, []).append(edge)
        edge_depth_keys[edge] = canonical_edge(*parse_path(key))
    return incoming, outgoing, edge_depth_keys
