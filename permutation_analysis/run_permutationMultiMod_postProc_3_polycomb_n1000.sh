#!/bin/bash
##SBATCH --time=06:00:00
module load R/3.5.2
R --no-save --no-restore <<EOS
rmarkdown::render(
	'permutationMultiMod_postProc_3.Rmd',
	output_file = '../out/permutationMultiMod_postProc_3_polycomb_n1000.html',
	params = list(
		inFile = "Projects/tRNA_Stuff/EpiTwin/new_permutation/out/out_gitIg/polycomb_CpGd_perm_t45_N1000.Rds"
	)
)
EOS

