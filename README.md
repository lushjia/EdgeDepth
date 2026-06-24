# EdgeDepth

Trait association analysis using variants (edges) from a human pangenome reference.
This repository contains a set of Nextflow pipelines that take short-read WGS data,
align it to a pangenome graph, and produce a filtered, normalized, allele balance replaced per-edge depth
matrix ready for downstream trait association analysis.

The pipeline is organized as five sequential steps, each in its own numbered folder.
Each step's output feeds directly into the next step's input — run them in order.

## Pipeline overview

```
1_alignment_count
    WGS CRAM align to pangenome graph (vg giraffe)
        -> raw edge depth x sample matrix

2_normalization
    size-factor normalization across samples
        -> normalized edge depth x sample matrix

3_redundancy_filter
    select representative edges (removes redundant edges)
        -> representative edge depth (norm) x sample matrix

4_variable_edge_filter
    alignment depth/GMM/outlier filters keep only variable edges
        -> variable edge depth (norm) x sample matrix

5_biallelic_AB
    for simple biallelic variants, replace edge depth with allele balance
        -> variable edge [depth (norm) or allele balance] x sample matrix
```


## Requirements
- [Nextflow](https://www.nextflow.io/) (DSL2)
- A Slurm cluster profile (`mccleary`) and container profiles (`docker`/`singularity`/`apptainer`) are
  defined in each step's `nextflow.config` — adjust `executor`/`queue`/`clusterOptions` for your own cluster.

## How to run
Each step folder has its own `README.md` with that step's specific run command and information.

## Output

The final output of the full pipeline is edge by sample matrix
edges pass all filters with normalzied edge depth or allele balance if the variants is a biallelic variant. 
convert to vcf 

## Toy example

A small toy dataset (2-sample CRAMs, a small reference region, and a matching pangenome
GBZ/HAPL/edges file) will be added to `data/` for quickly testing the pipeline end to end.
