# Hard genotype (optional)
This is an optional step that converts raw edge depth into discrete (hard) genotypes for every edge. It is run after
3_redundancy_filter and is independent of steps 4 and 5.

Each edge is genotyped with one of two models:
* **Allele (multinomial) model**: within each region of the graph bounded by two GRCh38 nodes, the region is turned
  into a directed acyclic graph. Starting with ploidy 2, the ploidy of every node is split across its outgoing edges
  by choosing the most likely genotype given the edges' read depths. A node's ploidy is the sum of the dosages of its
  incoming edges. The edge genotype is its dosage (0, 1, 2, ...).
* **Copy-number (Poisson) model**: for edges inside directed cycles, edges whose depth deviates from the allele call,
  and edges outside the regions above. The genotype is the most likely integer copy number given the edge depth and
  the local per-haplotype depth (lambda) estimated from nearby GRCh38 nodes.

When an edge has both calls, the allele dosage is reported.

The pipeline outputs:

* An edge-by-sample genotype matrix for each chromosome.
* Per-sample call details (genotype, depths, likelihoods, GQ, PL) for every call.

## Pipeline overview
This pipeline performs the following steps:

1. Per chromosome (sample-independent), using the GFA and the edge list:
    * Get the GRCh38 node order from the GRCh38 walk in the GFA.
    * Build the edge graph, remove bridge edges that are not on GRCh38, and split the graph into biconnected
      components at articulation nodes.
    * For components bounded by two GRCh38 nodes, build a directed graph from the haplotype walks in the GFA,
      collapse directed cycles to get a DAG, and record the cycle edges.
    * Components with only one GRCh38 boundary node are genotyped with the copy-number model.
    * Edges that cannot be genotyped this way (components with no GRCh38 boundary node, removed bridge edges, and
      unreachable DAG edges) are grouped by connected component. Each group gets the nearest GRCh38 node as its
      anchor for lambda and is genotyped with the copy-number model.
2. Per chromosome and per batch of samples, genotype every edge of each sample:
    * Estimate lambda for each component: the median depth at the nearest ~20 GRCh38 articulation nodes, divided by 2.
    * Run the allele model over each component's DAG.
    * Flag allele-called edges whose depth deviates from `lambda x dosage` (two-sided Poisson p < 0.01, unless
      more than half of the component's edges deviate), and run the copy-number model on them, on cycle edges,
      and on the edges from the last bullet of step 1.
    * Final genotype: the allele dosage if available (uncalled allele sites give 0), otherwise the copy number,
      or `.` if lambda could not be estimated.
3. Merge the per-sample genotypes into an edge-by-sample matrix per chromosome.

## How to run the pipeline
```
nextflow run hard_genotype.nf \
    -profile      mccleary \
    --chroms      chr21 \ # default chr1-22
    --gfa_dir     data/ \
    --edge_raw_dir ../3_redundancy_filter/results \
    --depth       ../2_normalization/results/all_sample.hprc-v2.0-mc-grch38.edge_depth.txt \
    --batch_size  10 \
    --scripts_dir scripts \
    --outdir      results
```
By default, `--chroms` runs on chromosomes 1–22. To run a subset of chromosomes, provide a single chromosome, such as chr21, or a comma-separated list, such as chr21,chr22.      
Add `-resume` to resume from cached work after a failed or interrupted run.    

## Required input data
| Parameter | Input type | Description |
|:----|:------|:----------|
| `--gfa_dir` | Directory | Directory containing the pangenome graph file in GFA format for each <br>chromosome. Files should be named using the format chr*.gfa. Each file <br>must contain exactly one GRCh38 walk (W-line). |
| `--edge_raw_dir` | Directory | Directory containing per-chromosome raw edge-depth files, named <br>chr*.depth_per_edge.txt.gz. These files are generated in <br>step 3_redundancy_filter. |
| `--depth` | File | Genome-wide raw edge-depth matrix generated in step 2_normalization <br>(the `--depth` input of step 3_redundancy_filter). Only its header line is <br>read, to get the sample names (columns 3 onward). |
| `--scripts_dir` | Scripts directory | Directory containing the Python scripts used in this step.|
| `--outdir` | Directory | Directory for output files. |

## Additional editable parameters
| Parameter | Input type | Description |
|:----|:------|:----------|
| `--chroms` | Chromosome list | Chromosomes to process. Provide a single chromosome, such as chr21, <br>or a comma-separated list, such as chr21,chr22. If not specified, <br>the default is chromosomes 1–22.|
| `--batch_size` | Integer | Number of samples genotyped per job (default 10). Each job reads the <br>chromosome's depth file once for the whole batch. Memory scales with <br>batch size (~30 GB per 10 samples on large chromosomes). |

More details about each step and parameter can be found in the Nextflow pipeline.

## Output files
The pipeline produces the following outputs:         
1. Edge genotype matrix (`chr*.edge_genotype_matrix.tsv`)

An edge-by-sample matrix (`edge_id`, one column per sample) with the hard genotype of every edge in every sample:
allele dosage or copy number. `.` means no call (no local depth estimate). Edge ids are in canonical orientation
(an edge and its reverse complement, e.g. `>13>14` and `<14<13`, are the same edge, reported starting with `>` when possible).

2. Per-sample call details (`genotype_info/chr*/chr*.{sample}.edge_genotype_info.tsv`)

One row per call:

| Column | Description |
|:----|:----------|
| `edges` | Allele call: comma-separated edges leaving the same node (the alleles, in GT index order). <br>CNV call: a single edge. |
| `CALL_TYPE` | `Allele` (multinomial model) or `CNV` (copy-number model). An edge can have both; <br>the matrix reports the allele dosage. |
| `FORMAT` | `GT:AD:DP:GL:GQ:PL` |
| `SAMPLE` | GT: allele genotype (e.g. `0/1`, `.` if uncalled) or copy number; AD: edge depth(s); <br>DP: total depth; GL: log10 likelihoods; GQ: genotype quality; PL: phred-scaled likelihoods. |
