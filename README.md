# EdgeDepth

Trait association analysis using variants (edges) from a human pangenome reference.
This repository contains a set of Nextflow pipelines that take short-read WGS data,
align it to a pangenome graph, and produce a filtered, normalized, per-edge depth
matrix suitable for downstream association testing.

The pipeline is organized as five sequential steps, each in its own numbered folder.
Each step's output feeds directly into the next step's input — run them in order.

## Pipeline overview

```
1_alignment_count      WGS CRAM -> align to pangenome graph (vg giraffe) -> per-sample edge depth
        |
2_normalization        merge per-sample depth into one matrix -> DESeq2-style size-factor normalization
        |
3_redundancy_filter     split matrix by chromosome -> bubble-puncture filter removes redundant edges
        |
4_variable_edge_filter  depth/GMM/4-MAD filters keep only genuinely variable edges
        |
5_biallelic_AB          for simple biallelic snarls, replace edge depth with allele balance
```

| Step | Folder | What it does | Key output |
|---|---|---|---|
| 1 | `1_alignment_count` | Align CRAM/FASTQ to the pangenome graph with `vg giraffe`, count read depth per edge | `{sample}.hprc-v2.0-mc-grch38.edge_depth.txt` (per sample) |
| 2 | `2_normalization` | Combine per-sample depth files into one matrix; compute and apply DESeq2-style median-of-ratios size factors | `all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt` (genome-wide) |
| 3 | `3_redundancy_filter` | Split depth matrices by chromosome; remove edges made redundant by graph topology (bubble-puncture method) | `{chrom}.independent_edge_list.normalized_depth.txt` |
| 4 | `4_variable_edge_filter` | Depth thresholds + GenomeSTRiP GMM + 4-MAD rescue to keep only variable (non-monomorphic) edges | `{chrom}.variable_edge_list.normalized_depth.txt` |
| 5 | `5_biallelic_AB` | For triangle/diamond snarls on the GRCh38 path, replace the kept edge's depth with its allele balance vs. the GRCh38 edge | `{chrom}.variable_edge_list.normalized_depth.AB_replaced.txt` |

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- `vg` (graph alignment, edge depth counting)
- `samtools` (CRAM/FASTQ conversion)
- `kmc` (k-mer counting for haplotype sampling)
- `R` (GenomeSTRiP-style GMM fitting in step 4)
- Python 3 with `pandas`, `numpy`, `scipy`, `networkx`, `pysam`
- A Slurm cluster profile (`mccleary`) and container profiles (`docker`/`singularity`/`apptainer`) are
  defined in each step's `nextflow.config` — adjust `executor`/`queue`/`clusterOptions` for your own cluster.

## Repository layout

```
1_alignment_count/      align_count_edgedepth.nf, nextflow.config
2_normalization/        normalize_edgedepth.nf, nextflow.config, scripts/
3_redundancy_filter/    redundancy_filter.nf, nextflow.config, script/
4_variable_edge_filter/ variable_edge_filter.nf, nextflow.config, scripts/
5_biallelic_AB/         biallelic_AB.nf, nextflow.config, scripts/
data/                   shared reference files (edges list, zjoin, toy example data — TBD)
```

Each step folder also has its own `README.md` with that step's specific run command and parameters.

## How to run

Run each step from inside its own folder, pointing `--outdir` of one step at the input directory
parameter of the next. All five pipelines are run per chromosome from step 3 onward.

```bash
ml Nextflow/25.04.6

# Step 1: align and count edge depth per sample
cd 1_alignment_count
nextflow run align_count_edgedepth.nf -profile mccleary \
    --cram_list data/samples_insertsize.txt \
    --b38_ref   /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --gbz       /path/to/hprc-v2.0-mc-grch38.gbz \
    --hapl      /path/to/hprc-v2.0-mc-grch38.hapl \
    --edges     ../data/hprc-v2.0-mc-grch38.edges.txt \
    --zjoin     ../data/zjoin \
    --outdir    results

# Step 2: merge + normalize
cd ../2_normalization
nextflow run normalize_edgedepth.nf -profile mccleary \
    --samples_list data/samples_list.txt \
    --edges        ../data/hprc-v2.0-mc-grch38.edges.txt \
    --depth_dir    ../1_alignment_count/results \
    --outdir       results

# Step 3: redundancy filter (bubble puncture), run chr1-chr22
cd ../3_redundancy_filter
nextflow run redundancy_filter.nf -profile mccleary \
    --edges      ../data/hprc-v2.0-mc-grch38.edges.txt \
    --depth      ../2_normalization/results/all_sample.hprc-v2.0-mc-grch38.edge_depth.txt \
    --norm_depth ../2_normalization/results/all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt \
    --gfa_dir    /path/to/per_chrom_gfa \
    --outdir     results

# Step 4: variable edge filter (GMM + 4-MAD)
cd ../4_variable_edge_filter
nextflow run variable_edge_filter.nf -profile mccleary \
    --norm_depth_dir ../3_redundancy_filter/results \
    --gmmscript      /path/to/compute_gmm.R \
    --outdir         results

# Step 5: replace biallelic edges with allele balance
cd ../5_biallelic_AB
nextflow run biallelic_AB.nf -profile mccleary \
    --vcf_dir        /path/to/per_chrom_vcf \
    --gfa_dir        /path/to/per_chrom_gfa \
    --snarls         /path/to/hprc-v2.0-mc-grch38.snarls.json \
    --edge_raw_dir   ../3_redundancy_filter/results \
    --kept_edges_dir ../3_redundancy_filter/results \
    --norm_depth_dir ../4_variable_edge_filter/results \
    --outdir         results
```

Add `-resume` to any of the above to resume from cached work after a failed or interrupted run.

## Output

The final output of the full pipeline is one file per chromosome:

```
{chrom}.variable_edge_list.normalized_depth.AB_replaced.txt
```

Tab-separated, no header: `chr  edge_id  sample1  sample2  ...  sampleN`. Each value is either the
normalized depth for that sample at that edge, or — for edges resolved by step 5 — the allele
balance of the kept edge against its GRCh38 alternative.

## Toy example

A small toy dataset (3-sample CRAMs, a small reference region, and a matching pangenome
GBZ/HAPL/edges file) will be added to `data/` for quickly testing the pipeline end to end.
