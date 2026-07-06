# Normalize edge depth across samples 
This pipeline combines the per-sample edge-depth files generate in step 1 and normalizes raw edge depth 
across samples to account for sequencing depth and library-size differences.

The pipeline outputs two edge-by-sample matrices: 
- one containing raw edge depth
- one containing normalized edge depth.

## Pipeline overview
This pipeline performs the following steps:

1. Combine per-sample edge depth files into a single raw edge-by-sample matrix.
2. Compute per-sample size factors using the DESeq2 median-of-ratios method in a streaming implementation.
3. Normalize edge depth across samples using the estimated size factors.

## How to run the pipeline 
```
nextflow run align_count_edgedepth.nf \
	-profile mccleary \
	--edges hprc-v2.0-mc-grch38.edges.txt \
	--depth_dir ../1_alignment_count/results \
	--scripts_dir scripts \
	--outdir results

```
Add `-resume`to resume from cached work after a failed or interrupted run.


## Required input data
| Parameter | Input type | Description |
|:----|:------|:----------|
| `--edges` | Edge list | Tab-separated file containing chromosome and edge IDs. <br>This file defines the edge order in the step 1 output. |
| `--depth_dir` | Directory | Directory containing the per-sample *.edge_depth.txt files generated in Step 1. |
| `--scripts_dir` | scripts directory | Directory containing the Python scripts used for normalization. |
| `--outdir`| Directory | Directory for output files. |

`--edges`: The edge list file is the same as the one used in Step 1.

More details about each step and parameter can be found in the Nextflow pipeline.

## Intermediate files 
Intermediate files are shown as symbolic links in the output directory. They are stored at work directory, useful for checking intermediate results, but are not used in subsequent steps and can be deleted.		

The size-factor file `all_sample.size_factors.txt` contains one line per sample, with the sample ID and the corresponding size factor. This file is used to normalize edge depth for each sample.

## Output file
The pipeline outputs two edge-by-sample matrices:
1. Raw edge-depth matrix (`all_sample.hprc-v2.0-mc-grch38.edge_depth.txt`): Matrix containing raw read depth for each edge in each sample.
2. Normalized edge-depth matrix (`all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt`): Matrix containing normalized edge depth after adjusting for per-sample sequencing depth and library-size differences.

In both matrices, rows represent graph edges and columns represent samples.
