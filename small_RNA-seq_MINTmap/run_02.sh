#!/bin/bash
##SBATCH --time=06:00:00
module load R/3.5.2
R --no-save --no-restore <<EOS
rmarkdown::render(
	'tRNAfragPermR_multiplan_1.Rmd',
	output_file = 'tRNAfragPermR_multiplan_1.html',
	params = list(
		pwd = "/scratch/rja1e16/PhD_git/Projects/tRNA_Stuff/tRNA_expression/src/",
		outDir = "/scratch/rja1e16/PhD_git/Projects/tRNA_Stuff/tRNA_expression/out_run/"
	)
)
EOS

