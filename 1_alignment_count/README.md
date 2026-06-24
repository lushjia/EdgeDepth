# WGS alignment to a pangenome reference
This pipeline can be used to perform short read alignmetn to the human pangenome reference with vg giraffe.

## What the pipeline can do
index graph, haplotype sampling, alignment 

## How to run the pipeline
nextflow run main.nf 

## Required input data
fastq 
gfa/gbz 

## output format 





# how to use

# script and config - copy to your folder: 
	main.nf
	nextflow.config

# input - in data folder - format: header(fixed) + 1 sample per row; 3 columns; tab seperator:
cram	insert_size	std
/gpfs/gibbs/pi/ycgh/ih233/projects/METSIM_120525/cram/2893844947/2893844947.cram	360.0	77.27152809609369

# run 
# compute node(2cpu, 5G, time=28-00:00:00)
# current node limit: 20

ml Nextflow/25.04.6 

nextflow run main.nf \
	-resume \
	-ansi-log true \
	-profile mccleary \
	--cram_list data/samples_insertsize.txt \
	--b38_ref /gpfs/gibbs/pi/ycgh/lushjia/reference/hg38_alt/GRCh38_full_analysis_set_plus_decoy_hla.fa \
	--gbz /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/graph/hprc-v2.0-mc-grch38.gbz \
	--hapl /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/graph/hprc-v2.0-mc-grch38.hapl \
	--outdir results

# output 
results will be in results folder
intermediate files will be in work folder - can be deleted later
output edge count for each sample
