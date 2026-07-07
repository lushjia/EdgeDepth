# EdgeDepth

EdgeDepth is a Nextflow pipeline for measuring genotype information of pangenome graph-based variants from short-read sequencing data using graph features (edges). The resulting edge depth matrix can be used in downstream trait association analyses, including eQTL, caQTL, and other molecular or complex trait studies (GWAS). 

**Input**: short-read WGS + pangenome graph        
**Output**: A filtered edge-by-sample matrix containing filtered, normalized edge-depth values, or allele-balance values for edges represented biallelic variants

## Pipeline overview

```
1_alignment_count
    align WGS to pangenome graph
        -> raw edge depth count x sample matrix

2_normalization
    normalize edge depth across samples
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

## Requirements
- [Nextflow](https://www.nextflow.io/) (DSL2)
- singularity
- A cluster profile (`mccleary`) — adjust `executor`/`queue`/`clusterOptions` for your own cluster in `nextflow.config`.

## How to run + toy example
Each step folder has its own `README.md` with that step's specific run command and information.

A small toy dataset is provided in the `data/` directory for quickly testing the pipeline end to end.

Use the following commands to run the toy example:

```bash
git clone https://github.com/lushjia/EdgeDepth.git

cd EdgeDepth

# step 1 
nextflow run 1_alignment_count/align_count_edgedepth.nf \
  -resume \
  -profile mccleary \
  --cram_list data/samples_insertsize.txt \
  --b38_ref data/chr21.GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --gbz data/toy.gbz \
  --hapl data/toy.hapl \
  --edges data/toy.edges.tsv \
  --scripts_dir 1_alignment_count/scripts \
  --outdir results

# step 2  
nextflow run 2_normalization/normalize_edgedepth.nf \
    -resume \
    -profile mccleary \
    --edges        data/toy.edges.tsv \
    --depth_dir    results \
    --scripts_dir  2_normalization/scripts \
    --outdir       results

# step 3
nextflow run 3_redundancy_filter/redundancy_filter.nf \
    -resume \
    -profile mccleary \
    --chroms    chr21 \
    --edges      data/toy.edges.tsv \
    --depth      results/all_sample.hprc-v2.0-mc-grch38.edge_depth.txt \
    --norm_depth results/all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt \
    --gfa_dir    data \
    --scripts_dir 3_redundancy_filter/scripts \
    --outdir     results

# step 4 
nextflow run 4_variable_edge_filter/variable_edge_filter.nf \
    -resume \
    -profile mccleary \
    --norm_depth_dir results \
    --scripts_dir    4_variable_edge_filter/scripts \
    --outdir         results \
    --min_samples    0 \
    --min_depth      0 \
    --gmm_min_outlier_samples 0\
    --mad_min_outlier_samples 0

# step 5
nextflow run 5_biallelic_AB/biallelic_AB.nf \
    -resume \
    -profile         mccleary \
    --chroms         chr21 \
    --vcf_dir        data \
    --gfa_dir        data \
    --snarls         data/toy.snarls.json \
    --edge_raw_dir   results \
    --kept_edges_dir results \
    --norm_depth_dir results \
    --scripts_dir    5_biallelic_AB/scripts \
    --outdir         results
```




## Citation: 
If you use `EdgeDepth` pipeline in your work, please cite:
> S. Lu, W.-W. Liao, M. K. DeGorter, P. C. Goddard, J. Ebler, T.-Y. Lu, Human Pangenome Reference Consortium, M. J. P. Chaisson, T. Marschall, S. B. Montgomery, N. O. Stitziel, and I. M. Hall.         
> Pangenome-based human genome analysis improves trait association and genomic prediction        
> bioRvix,
> doi: https://doi.org/10.64898/2026.07.01.735728

