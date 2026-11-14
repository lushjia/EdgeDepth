"""
Genotyping models shared by the 6_hard_genotype_optional scripts.

1. Allele (multinomial) model - for edges leaving one node in a component's
   DAG: the node's ploidy is split across its outgoing edges; every candidate
   genotype is scored with a multinomial likelihood of the edge read depths
   (symmetric error rate, svtyper-style GQ/PL).
2. Copy-number (Poisson) model - for cyclic edges, edges whose depth deviates
   from the allele call, and recovered edges: integer copy number scored
   against Poisson(lambda * cn), lambda = local per-haplotype depth.
3. Local depth (lambda) - median depth at the ~20 GRCh38 articulation nodes
   nearest the component, divided by 2 (diploid).
"""

import bisect
import itertools
import math
from collections import deque
from dataclasses import dataclass
from typing import List, Optional, Tuple

from scipy.stats import poisson

from graph_utils import canonical_edge, parse_path

Genotype = Tuple[int, ...]   # sorted allele indices, e.g. (0, 1) for "0/1"
_LN10 = math.log(10.0)


# ---------------------------------------------------------------------------
# 1. Allele (multinomial) model
# ---------------------------------------------------------------------------
@dataclass
class GenotypeCall:
    """Result of genotyping one site (the outgoing edges of one node)."""
    genotype: Optional[Genotype]          # None when not callable
    gq: float                             # phred gap between best and runner-up
    genotypes: List[Genotype]             # all candidates considered
    log10_likelihoods: List[float]        # parallel to `genotypes`
    log10_posteriors: List[float]         # normalized, flat prior
    pl: List[int]                         # phred-scaled likelihoods, best genotype = 0
    depth: float                          # total observations at the site

    def as_vcf_string(self):
        """VCF-style GT, e.g. '0/1'; '.' when uncalled."""
        if self.genotype is None:
            return "."
        return "/".join(str(allele) for allele in self.genotype)


def symmetric_fractions(genotype, n_alleles, error_rate=0.001):
    """
    Expected fraction of reads per allele given a genotype; alleles absent from
    the genotype share `error_rate` equally.
    e.g. (1, 1), 3 alleles -> [0.001/2, 1 * (1 - 0.001), 0.001/2]
    """
    ploidy = len(genotype)
    counts = [0] * n_alleles
    for allele in genotype:
        counts[allele] += 1
    n_absent = sum(1 for c in counts if c == 0)
    if n_absent == 0:
        return [c / ploidy for c in counts]
    per_absent = error_rate / n_absent
    return [per_absent if c == 0 else (c / ploidy) * (1.0 - error_rate) for c in counts]


def log10_multinomial_coefficient(counts):
    """log10( n! / (a_1! ... a_k!) ), via lgamma (works for non-integer counts)."""
    total = math.fsum(counts)
    log_coeff = math.lgamma(total + 1.0) - math.fsum(math.lgamma(c + 1.0) for c in counts)
    return log_coeff / _LN10


def multinomial_log10_pmf(counts, probs, include_coefficient=True):
    """log10 probability of `counts` given category probabilities."""
    total = 0.0
    for count, prob in zip(counts, probs):
        if count == 0:
            continue
        total += count * math.log10(prob)
    if include_coefficient:
        total += log10_multinomial_coefficient(counts)
    return total


def _log10_sum_exp(values):
    """Numerically stable log10(sum(10**v))."""
    finite = [v for v in values if v != float("-inf")]
    if not finite:
        return float("-inf")
    largest = max(finite)
    return largest + math.log10(math.fsum(10.0 ** (v - largest) for v in finite))


def bayes_gt(counts, ploidy=2, error_rate=0.001, include_coefficient=True):
    """
    Log10 likelihood of every genotype (VCF GT/GL/PL order) given allele counts.
    e.g. bayes_gt([5, 5]) -> ([(0, 0), (0, 1), (1, 1)], [-12.60, -0.61, -12.60])
    """
    n_alleles = len(counts)
    genotypes = sorted(
        itertools.combinations_with_replacement(range(n_alleles), ploidy),
        key=lambda gt: gt[::-1],
    )
    likelihoods = [
        multinomial_log10_pmf(counts, symmetric_fractions(gt, n_alleles, error_rate), include_coefficient)
        for gt in genotypes
    ]
    return genotypes, likelihoods


def call_genotype(counts, ploidy=2, error_rate=0.001, min_depth=1.0, max_gq=200.0, max_pl=255):
    """
    Maximum-likelihood genotype with GQ (phred gap to the 2nd best, capped at
    max_gq) and PL. Uncalled when ploidy is 0 or total depth < min_depth.
    """
    depth = math.fsum(counts)
    if ploidy == 0 or depth < min_depth:
        return GenotypeCall(None, 0.0, [], [], [], [], depth)

    genotypes, likelihoods = bayes_gt(counts, ploidy, error_rate, include_coefficient=True)

    normalizer = _log10_sum_exp(likelihoods)
    posteriors = (
        [lk - normalizer for lk in likelihoods]
        if normalizer != float("-inf")
        else [float("-inf")] * len(likelihoods)
    )
    best_lk = max(likelihoods)
    pl = [min(max_pl, max(0, round(-10.0 * (lk - best_lk)))) for lk in likelihoods]

    if len(genotypes) < 2:
        return GenotypeCall(genotypes[0], 0.0, genotypes, likelihoods, posteriors, pl, depth)

    ranked = sorted(range(len(likelihoods)), key=lambda i: likelihoods[i], reverse=True)
    best, second = ranked[0], ranked[1]
    gq = min(-10.0 * (likelihoods[second] - likelihoods[best]), max_gq)
    return GenotypeCall(genotypes[best], gq, genotypes, likelihoods, posteriors, pl, depth)


def genotype_function(edge_counts, node_ploidy):
    """
    Genotype the outgoing edges of one node.
    edge_counts : {edge: depth}
    Returns ({edge: dosage}, {tuple(edges): GenotypeCall}). Uncalled -> dosage 0
    for every edge (which also passes ploidy 0 downstream).
    """
    edge_ids = list(edge_counts.keys())
    counts = [edge_counts[e] for e in edge_ids]
    call = call_genotype(counts, node_ploidy)

    dosage = [0] * len(edge_ids)
    if call.genotype is not None:
        for allele_index in call.genotype:
            dosage[allele_index] += 1
    return dict(zip(edge_ids, dosage)), {tuple(edge_ids): call}


def genotype_dag(incoming, outgoing, edge_depth_keys, seed_node, seed_ploidy, edge_depths):
    """
    Breadth-first flow genotyping over one component's DAG. Starting from the
    seed node (ploidy 2), the outgoing edges of each node are genotyped together
    with the node's ploidy; a node's ploidy is the sum of the dosages of its
    incoming edges, available once all of them are genotyped.
    """
    node_ploidy = {seed_node: seed_ploidy}
    edge_genotype = {}
    edgegroup_info = {}
    queue = deque([str(seed_node)])
    processed_nodes = set()
    while queue:
        node = queue.popleft()
        if node in processed_nodes:
            continue
        node_outgoing = outgoing.get(node, [])
        untyped_edges = [e for e in node_outgoing if e not in edge_genotype]
        if untyped_edges:
            depths = {e: edge_depths.get(edge_depth_keys[e], 0) for e in untyped_edges}
            new_geno, new_info = genotype_function(depths, node_ploidy[node])
            edge_genotype.update(new_geno)
            edgegroup_info.update(new_info)
        processed_nodes.add(node)
        for edge in node_outgoing:
            target = edge[1]
            target_incoming = incoming.get(target, [])
            if all(e in edge_genotype for e in target_incoming) and target not in node_ploidy:
                node_ploidy[target] = sum(edge_genotype[e] for e in target_incoming)
                queue.append(target)
    return edge_genotype, edgegroup_info


# ---------------------------------------------------------------------------
# 2. Copy-number (Poisson) model
# ---------------------------------------------------------------------------
def call_copy_number(observed_depth, lam, max_cn=None, error_rate=0.001, max_gq=200.0, max_pl=255):
    """
    Maximum-likelihood copy number for one edge under Poisson(lam * cn)
    (cn = 0 uses lam * error_rate). max_cn defaults to depth/lam + 4 (min 4).
    Returns (best_cn, gq, copy_numbers, log10_likelihoods, pl).
    """
    if max_cn is None:
        max_cn = max(4, math.ceil(observed_depth / lam) + 4) if lam > 0 else 4

    copy_numbers = list(range(0, max_cn + 1))
    log10_likelihoods = [
        poisson.logpmf(observed_depth, lam * cn if cn > 0 else lam * error_rate) / _LN10
        for cn in copy_numbers
    ]
    ranked = sorted(range(len(copy_numbers)), key=lambda i: log10_likelihoods[i], reverse=True)
    best, second = ranked[0], ranked[1]
    best_ll = log10_likelihoods[best]

    if math.isfinite(best_ll):
        gq = min(max_gq, -10.0 * (log10_likelihoods[second] - best_ll))
        pl = [min(max_pl, max(0, round(-10.0 * (ll - best_ll)))) for ll in log10_likelihoods]
    else:
        # lam == 0 with observed depth > 0: every copy number impossible
        gq = 0.0
        pl = [0] * len(log10_likelihoods)
    return copy_numbers[best], gq, copy_numbers, log10_likelihoods, pl


def copy_number_info(observed_depth, lam):
    """call_copy_number packed as the dict written to the info file."""
    best_cn, gq, copy_numbers, log10_likelihoods, pl = call_copy_number(observed_depth, lam)
    return {
        "cn": best_cn,
        "depth": observed_depth,
        "lambda": lam,
        "copy_numbers": copy_numbers,
        "log10_likelihoods": log10_likelihoods,
        "gq": gq,
        "pl": pl,
    }


def genotype_component(entry, stem_depth, edge_depths, deviance_p=0.01, max_deviant_fraction=0.5):
    """
    Genotype one component for one sample.

      1. Allele model over the component's DAG (if it has one).
      2. Allele-called edges whose depth deviates from stem_depth * dosage
         (two-sided Poisson p < deviance_p) are flagged as possible CNV, unless
         more than max_deviant_fraction of the component's edges deviate (then
         lambda itself is assumed off and nothing is flagged).
      3. Copy-number model for cyclic edges + flagged edges (skipped if
         stem_depth is None).

    Returns:
        edge_genotype  : {canonical edge: dosage} from the allele model
        edgegroup_info : {tuple(raw edge keys): GenotypeCall}
        cnv_info       : {canonical edge: copy_number_info dict}
    """
    if entry["incoming"] is not None:
        dag_genotype, edgegroup_info = genotype_dag(
            entry["incoming"], entry["outgoing"], entry["edge_depth_keys"],
            entry["dag_seed_node"], entry["start_component_node_ploidy"], edge_depths,
        )
    else:
        dag_genotype, edgegroup_info = {}, {}

    cnv_edges = set(entry["cyclic_edges"])
    if stem_depth is not None:
        deviant_edges = set()
        for edge, dosage in dag_genotype.items():
            observed = edge_depths.get(canonical_edge(*parse_path(edge[2])), 0)
            expected = stem_depth * dosage if dosage > 0 else stem_depth * 0.001
            if observed <= expected:
                p_value = poisson.cdf(observed, expected)
            else:
                p_value = poisson.sf(observed - 1, expected)
            if p_value < deviance_p:
                deviant_edges.add(edge)
        if deviant_edges and len(deviant_edges) <= max_deviant_fraction * len(dag_genotype):
            cnv_edges.update(canonical_edge(*parse_path(e[2])) for e in deviant_edges)

    cnv_info = {}
    if stem_depth is not None:
        for edge in cnv_edges:
            observed = edge_depths.get(canonical_edge(*parse_path(edge)), 0)
            cnv_info[edge] = copy_number_info(observed, stem_depth)

    edge_genotype = {canonical_edge(*parse_path(e[2])): d for e, d in dag_genotype.items()}
    return edge_genotype, edgegroup_info, cnv_info


# ---------------------------------------------------------------------------
# 3. Local depth (lambda)
# ---------------------------------------------------------------------------
def compute_stem_total_depth(articulation_nodes, outgoing_paths, incoming_paths, end_node, edge_depths):
    """
    Depth at each articulation node = sum of depths of the edges leaving it
    (entering it, for the backbone end node). Only edges where the node is the
    '>'-oriented source (index 0) / target (index -1) are counted.
    """
    stem_total_depth = {}
    for node in articulation_nodes:
        if node == end_node:
            edges = {e for e in incoming_paths[node] if parse_path(e)[-1] == f">{node}"}
        else:
            edges = {e for e in outgoing_paths[node] if parse_path(e)[0] == f">{node}"}
        stem_total_depth[node] = sum(edge_depths.get(e, 0) for e in edges)
    return stem_total_depth


def local_lambda(stem_positions, path_edges, b38_node_order, start_component_node, end_component_node,
                 stem_total_depth, N=10, ploidy=2):
    """
    Per-haplotype depth around a component: median articulation-node depth over
    the component's articulation nodes plus N on either side, divided by ploidy.
    None if no articulation node is in the window.
    """
    i = bisect.bisect_left(stem_positions, b38_node_order[start_component_node])
    j = bisect.bisect_right(stem_positions, b38_node_order[end_component_node])
    lo = max(0, i - N + 1)
    hi = min(len(stem_positions), j + N)
    depths = sorted(stem_total_depth[path_edges[p]] for p in stem_positions[lo:hi])
    if not depths:
        return None
    return depths[len(depths) // 2] / ploidy
