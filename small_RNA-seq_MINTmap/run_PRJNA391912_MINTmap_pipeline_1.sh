#!/bin/bash
module load R/3.5.2
R --no-save --no-restore <<EOS
rmarkdown::render(
	'PRJNA391912_MINTmap_pipeline_1.Rmd',
	output_file = 'PRJNA391912_MINTmap_pipeline_1.html',
	params = list(
		fetchData = FALSE,
		plan = "plan(list(tweak(batchtools_slurm,resources = list(walltime = '00:15:00',ntasks=3L))))"
	)
)
EOS

