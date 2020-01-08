#!/bin/bash
##SBATCH --time=06:00:00
module load R/3.5.2
R --no-save --no-restore <<EOS
rmarkdown::render(
	'permutationMultiMod_3.Rmd',
	output_file = '../out/permutationMultiMod_3_bivalent_PromP_n1000.html',
	params = list(
		phdPath =  "/scratch/rja1e16/PhD_git/",
		nPerm = 1000L,
		features = "Projects/tRNA_Stuff/EpiTwin/data/wgEncodeAwgSegmentationChromhmmGm12878_PromP_ids_5col.bed",
		out = "Projects/tRNA_Stuff/EpiTwin/new_permutation/out/out_gitIg/bivalent-PromP_CpGd_perm_t45_N1000.Rds",
		modelData = 'Projects/tRNA_Stuff/EpiTwin/data/data_gitIg/BloodBatchCpGd3001_blacklistFiltered.Rds',
		threshold = 45L,
		plan = "plan(
			list(
				tweak(
					batchtools_slurm,
					resources = list(
						ntasks = 32L,
						walltime = '10:00:00',
						mem_per_cpu = 5000L
					)
				),
				tweak(
					multiprocess,
					workers = 32L
				)
			)
		)"
	)
)
EOS

