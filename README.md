# EdgeDepth

EdgeDepth is a pipeline for measuring genotype information of pangenome graph-based variants from short-read sequencing data using graph features (edges). The resulting edge depth matrix can be used in downstream trait association analyses, including eQTL, caQTL, and other molecular or complex trait studies (GWAS). 

**Input**: short-read WGS + pangenome graph        
**Output**: A filtered edge-by-sample matrix containing filtered, normalized edge-depth values, or allele-balance values for edges represented biallelic variants

## Pipeline overview

```
1_alignment_count
    WGS CRAM align to pangenome graph 
        -> raw edge depth x sample matrix

2_normalization
    normalized edge depth across samples
        -> normalized edge depth x sample matrix

3_redundancy_filter
    select representative edges (remove redundant edges)
        -> representative edge depth (norm) x sample matrix

4_variable_edge_filter
    keep only variable edges (alignment depth/GMM/outlier filters )
        -> variable edge depth (norm) x sample matrix

5_biallelic_AB
    replace edge depth with allele balance for simple biallelic variants
        -> variable edge [depth (norm) or allele balance] x sample matrix
```

To do: diagram 

## Requirements
- [Nextflow](https://www.nextflow.io/) (DSL2)
- singularity
- A cluster profile (`mccleary`) — adjust `executor`/`queue`/`clusterOptions` for your own cluster in `nextflow.config`.

## How to run
Each step folder has its own `README.md` with that step's specific run command and information.

## Toy example

A small toy dataset (2-sample CRAMs, a small reference region, and a matching pangenome
GBZ/HAPL/edges file) will be added to `data/` for quickly testing the pipeline end to end.
