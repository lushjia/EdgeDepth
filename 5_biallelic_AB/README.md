# Replace biallelic edge depth by allele balance
This pipeline replaces the normalized edge depth of biallelic edges with allele balance (AB). It is run after 4_variable_edge_filter. 
For biallelic edges, where edges represent biallelic variants, allele balance is calculated as:        
    
`AB = alternative allele edge raw depth / (alternative allele edge raw depth + reference allele edge raw depth)`        
    
The normalized edge-depth values for biallelic edges in the Step 4_variable_edge_filter edge-by-sample matrix are then 
replaced with the corresponding allele-balance values.

The pipeline outputs:

* An edge-by-sample matrix containing only biallelic edges, with allele balance values for each edge in each sample.    
* An edge-by-sample matrix containing the same edges as the Step 4_variable_edge_filter output,
  with allele balance for biallelic edges and normalized edge depth for multiallelic edges.

## Pipeline overview
This pipeline performs the following steps:

1. Map snarls to allele traversals/edges using the graph variant VCF.
2. For snarls without traversal information from VCF, recover recover allele traversals from GFA paths.
3. Traverse the snarl tree using parent-child relationships and retain only the edges first introduced at each snarl level, excluding edges already assigned to deeper child snarls.
4. Identify biallelic edges based on the edges present at each snarl level and the snarl structure.
Compute allele balance for each biallelic edge relative to the corresponding reference edge.
5. Replace normalized edge-depth rows for biallelic edges in the Step 4_variable_edge_filter variable-edge matrix with allele-balance rows. 

## How to run the pipeline
```
nextflow run biallelic_AB.nf -profile mccleary \
    --chroms         chr21 \ # default chr1-22 
    --vcf_dir        data/ \ 
    --gfa_dir        data/ \ # gfa file name: chr*.gfa
    --snarls         data/hprc-v2.0-mc-grch38.snarls.json \
    --edge_raw_dir   ../3_redundancy_filter/results \
    --kept_edges_dir ../3_redundancy_filter/results \
    --norm_depth_dir ../4_variable_edge_filter/results \
    --outdir         results
```
By default, `--chroms` runs on chromosomes 1–22. To run a subset of chromosomes, provide a single chromosome, such as chr21, or a comma-separated list, such as chr21,chr22.      
Add `-resume` to resume from cached work after a failed or interrupted run.    

## Required input data
| Parameter | Input type | Description |
|:----|:------|:----------|
| `--vcf_dir` | Directory | Directory containing graph variant VCF files that record the reference <br> and alternative traversal for each snarl. Files should be named using the <br>format hprc-v2.0-mc-grch38.raw.chr*.sorted.vcf.gz |
| `--gfa_dir` | Directory | Directory containing the pangenome graph file in GFA format for each <br>chromosome. Files should be named using the format chr*.gfa. |
| `--snarls` | JSON file | Snarl file recording the start and end nodes of each snarl, as well as <br>the start and end nodes of its parent snarl. |
| `--edge_raw_dir` | Directory | Directory containing per-chromosome raw edge-depth files, named <br>{chrom}.depth_per_edge.txt.gz. These files are generate in <br> step 3_redundancy_filter|
| `--kept_edges_dir` | Directory | Directory containing per-chromosome independent edge lists, named <br>{chrom}.independent_edge_list.txt. These files are generate <br>in step 3_redundancy_filter|
| `--norm_depth_dir` | Directory | Directory containing filtered edge-by-sample matrices, named <br>{chrom}.variable_edge_list.normalized_depth.txt. These files are <br>generated in step 4_biallelic_AB|
| `--outdir` | Directory | Directory for output files. |


## Additional editable parameters
| Parameter | Input type | Description |
|:----|:------|:----------|
| `--chroms` | Chromosome list | Chromosomes to process. Provide a single chromosome, such as chr21, <br>or a comma-separated list, such as chr21,chr22. If not specified, <br>the default is chromosomes 1–22.|

More details about each step and parameter can be found in the Nextflow pipeline.

The graph variant VCF files can be downloaded from HPRC: GRCh38_MinigraphCactus_variants CHM13_MinigraphCactus_variants    
The snarl file can be downloaded from HPRC: GRCh38_MinigraphCactus_snarl CHM13_MinigraphCactus_snarl    

gfa file can be converted from: link

## Output files
The pipeline produces the following outputs:         
1. Biallelic-edge allele-balance matrix

An edge-by-sample matrix  containing only biallelic edges, with allele balance values for each edge in each sample.

2. Final edge-by-sample matrix
An edge-by-sample matrix containing all variable edges retained after step 4_variable_edge_filter. In this matrix,
biallelic edges are represented by allele balance, while multiallelic edges are represented by normalized edge depth.
This file could be used for downstream trait association analysis. 




