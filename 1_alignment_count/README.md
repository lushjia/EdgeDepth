# Align WGS reads to a pangenome reference graph
This pipeline aligns short-read whole-genome sequencing (WGS) data in CRAM format to 
a human pangenome reference graph using vg giraffe with haplotype sampling. 
It then counts read depth for each graph edge and outputs one edge-depth file per sample.

Each output file contains one column of edge-depth values, with one row per graph edge. 
Edges are ordered according to the provided edge list, ensuring that output files 
from different samples can be directly combined into an edge-by-sample matrix.

## Pipeline overview 
This pipeline performs the following steps:  

1. Convert CRAM files to FASTQ using the provided reference FASTA.  
2. Align FASTQ reads to the pangenome graph (.gbz) using vg giraffe and the haplotype index file (.hapl).    
3. Count read depth for each graph edge and output one edge-depth file per sample, with edges ordered according to the provided edge list.  

## How to run the pipeline
```
nextflow run align_count_edgedepth.nf \
	-profile mccleary \
	--cram_list samples_insertsize.txt \
	--b38_ref GRCh38_full_analysis_set_plus_decoy_hla.fa \
	--gbz toy.gbz \
	--hapl toy.hapl \
	--edges toy.edges.tsv \
	--scripts_dir scripts \
	--outdir results_toy
```
Add `-resume`to resume from cached work after a failed or interrupted run.

## Required input data
| Parameter | Input type | Description |
|:----|:------|:----------|
| `-profile`| Configuration settings | Preset defined in the profiles block of nextflow.config that specifies configuration <br> settings for the execution environment. See: <br>[https://docs.seqera.io/nextflow/config#config-profiles](https://docs.seqera.io/nextflow/config#config-profiles) |
| `--cram_list` | Sample table | Path to a table listing short-read WGS CRAM files, <br>insert-size and corresponding standard deviation (std) information. |
| `--b38_ref` | FASTA | GRCh38 reference FASTA used to convert CRAM files to FASTQ. |
| `--gbz` | GBZ graph | Pangenome reference graph used for read alignment. |
| `--hapl` | Haplotype index | Haplotype index file used by `vg giraffe` for haplotype-aware read alignment. |
| `--edges` | Edge list | Tab-separated file containing chromosome and edge IDs. <br>This file defines the edge order in the output. | 
| `--scripts_dir` | Scripts directory | Directory containing the Python script (zjion) used to join edge-depth output <br>with the provided edge list, ensuring consistent edge order across samples. | 
| `--outdir`| Directory | Directory for output files. |

### Additional input details
- Example sample table
```
cram                      insert_size     std
cram/HG02976.final.cram   439             102.299400
cram/HG03130.final.cram   432             99.334200
```
- Example edge list
```
Chrom   Edge
chr10   >1>2
chr10   >2>3
chr10   >3>4
```

- The GRCh38 reference FASTA and associated index files can be downloaded from the 1000 Genomes Project: [link](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)

- The GBZ graph and HAPL haplotype index for the HPRC2 GRCh38 Minigraph-Cactus pangenome graph are available here: [`gbz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.gbz), [`hapl`](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.hapl)

- The edge list can be generated from the pangenome graph GFA file of each chromosome by the following command:
  
```
chr=chr1
awk -F'\t' -v chr=$chr '
BEGIN { OFS="\t"; print "Chrom","Edge" }
$1=="L" {
    node1=$2; strand1=$3; node2=$4; strand2=$5
    direct1 = (strand1=="+") ? ">" : "<"
    direct2 = (strand2=="+") ? ">" : "<"
    if (node1+0 < node2+0) {
        e1 = direct1 node1
        e2 = direct2 node2
    } else {
        rev_direct2 = (direct2==">") ? "<" : ">"
        rev_direct1 = (direct1==">") ? "<" : ">"
        e1 = rev_direct2 node2
        e2 = rev_direct1 node1
    }
    print chr, e1 e2
}
' $chr.gfa >> hprc-v2.0-mc-grch38.edges.tsv
```
  the pangenome graph VG file of each chromosome (GRCh38) can be found at: [pangenome graph per chr](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.chroms/);	            
  the pangenome graph GFA file can be converted from VG by [vg toolkit](https://github.com/vgteam/vg): `vg convert -f chr1.vg > chr1.gfa`

More details about each step and parameter can be found in the Nextflow pipeline.

## Output format
The pipeline outputs one file per sample. Each file contains a single column of read-depth values, 
with one row per graph edge. The row order follows the input edge list, allowing edge-depth files 
from multiple samples to be combined into a consistent edge-by-sample matrix in next step.

